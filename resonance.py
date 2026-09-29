"""Adaptive spectral feature search and optional Fano fit for finite FMM stacks."""
from dataclasses import asdict
from run_jobs import progress

import numpy as np
from scipy.optimize import curve_fit
from stack import StackModel, solve_stack, stack_convergence


def _fano(x, offset, slope, amplitude, x0, gamma, q):
    eps = 2*(x-x0)/gamma
    return offset + slope*(x-x0) + amplitude*(q+eps)**2/(1+eps**2)


def _lorentzian(x, offset, slope, amplitude, x0, gamma):
    return offset + slope*(x-x0) + amplitude/(1+(2*(x-x0)/gamma)**2)


def adaptive_resonance(model: StackModel, start: float, stop: float,
                       quantity: str = "R", extremum: str = "max",
                       points: int = 31, rounds: int = 3) -> dict:
    if quantity not in ("R", "T", "A") or extremum not in ("max", "min"):
        raise ValueError("Choose R, T, or A and a maximum or minimum")
    if not np.isfinite([start, stop]).all() or not 0 < start < stop:
        raise ValueError("Use positive increasing wavelength limits")
    if not 11 <= points <= 101 or not 1 <= rounds <= 5:
        raise ValueError("Use 11–101 points and 1–5 refinement rounds")
    samples = {}
    lo, hi = float(start), float(stop)
    fit_lo, fit_hi = lo, hi
    for round_index in range(rounds):
        progress(round_index, rounds, "Adaptive resonance rounds")
        fit_lo, fit_hi = lo, hi
        xs = np.linspace(lo, hi, points)
        local = []
        for wavelength in xs:
            current = StackModel(**{**asdict(model), "layers": model.layers,
                                    "wavelength_um": float(wavelength)})
            try:
                result = solve_stack(current)
                row = {"wavelength_um": float(wavelength),
                       **{key: float(result[key]) for key in ("R", "T", "A")}}
                samples[round(float(wavelength), 14)] = row
                local.append(row)
            except (ValueError, np.linalg.LinAlgError):
                continue
        if len(local) < 3:
            raise ValueError("Too few valid wavelengths for adaptive refinement")
        values = np.asarray([row[quantity] for row in local])
        i = int(np.argmax(values) if extremum == "max" else np.argmin(values))
        left, right = max(0, i-2), min(len(local)-1, i+2)
        if left == right:
            break
        lo, hi = local[left]["wavelength_um"], local[right]["wavelength_um"]
    rows = sorted(samples.values(), key=lambda row: row["wavelength_um"])
    vals = np.asarray([row[quantity] for row in rows])
    i = int(np.argmax(vals) if extremum == "max" else np.argmin(vals))
    feature = rows[i]
    refined = StackModel(**{**asdict(model), "layers": model.layers,
                            "wavelength_um": feature["wavelength_um"]})
    convergence = stack_convergence(refined)
    fit = None
    fit_error = None
    # Fit only the final refined window; covariance supplies approximate 1-sigma errors.
    fit_rows = [row for row in rows if fit_lo-1e-15 <= row["wavelength_um"] <= fit_hi+1e-15]
    if len(fit_rows) >= 11:
        x = np.asarray([row["wavelength_um"] for row in fit_rows])
        y = np.asarray([row[quantity] for row in fit_rows])
        span = x[-1]-x[0]
        try:
            popt, pcov = curve_fit(_fano, x, y,
                p0=[float(np.median(y)), 0, float((y.max()-y.min()) or .01),
                    feature["wavelength_um"], span/8, 1],
                bounds=([-2, -1e5, -5, x[0], span/10000, -100],
                        [2, 1e5, 5, x[-1], span*2, 100]), maxfev=30000)
            predicted = _fano(x, *popt)
            sigma = np.sqrt(np.maximum(0, np.diag(pcov)))
            fit = {"model": "Fano", "offset": float(popt[0]), "slope": float(popt[1]),
                   "amplitude": float(popt[2]), "lambda0_um": float(popt[3]),
                   "linewidth_um": float(abs(popt[4])), "fano_q": float(popt[5]),
                   "Q": float(abs(popt[3]/popt[4])),
                   "rmse": float(np.sqrt(np.mean((predicted-y)**2))),
                   "lambda0_sigma_um": float(sigma[3]), "linewidth_sigma_um": float(sigma[4]),
                   "wavelength_um": x.tolist(), "predicted": predicted.tolist()}
        except (RuntimeError, ValueError, FloatingPointError):
            # Symmetric resonances can make the Fano asymmetry parameter poorly
            # conditioned. A Lorentzian is the q→∞ limiting line shape.
            try:
                edge = float((y[0]+y[-1])/2)
                popt, pcov = curve_fit(_lorentzian, x, y,
                    p0=[edge, 0, float((y.max() if extremum == "max" else y.min())-edge),
                        feature["wavelength_um"], span/8],
                    bounds=([-2, -1e5, -5, x[0], span/10000],
                            [2, 1e5, 5, x[-1], span*2]), maxfev=30000)
                predicted = _lorentzian(x, *popt)
                sigma = np.sqrt(np.maximum(0, np.diag(pcov)))
                fit = {"model": "Lorentzian", "offset": float(popt[0]),
                       "slope": float(popt[1]), "amplitude": float(popt[2]),
                       "lambda0_um": float(popt[3]), "linewidth_um": float(abs(popt[4])),
                       "fano_q": None, "Q": float(abs(popt[3]/popt[4])),
                       "rmse": float(np.sqrt(np.mean((predicted-y)**2))),
                       "lambda0_sigma_um": float(sigma[3]), "linewidth_sigma_um": float(sigma[4]),
                       "wavelength_um": x.tolist(), "predicted": predicted.tolist()}
            except (RuntimeError, ValueError, FloatingPointError) as exc:
                fit_error = str(exc)
    if fit is not None:
        width, center = fit['linewidth_um'], fit['lambda0_um']
        spacing = float(np.min(np.diff(x)))
        edge = min(center-x[0], x[-1]-center)
        if (edge <= max(span*1e-4, width/2) or width >= span*1.999
                or width < 2*spacing or not np.isfinite(fit['linewidth_sigma_um'])
                or fit['linewidth_sigma_um'] >= width):
            fit_error = ('Unresolved fit: the center/linewidth is constrained by the window, '
                         'undersampled, or uncertain. Widen the interval to include both flanks '
                         'and background, increase spectral sampling, and refit. Q is withheld.')
            fit = None
    return {"rows": rows, "quantity": quantity, "extremum": extremum,
            "feature": feature, "fit": fit, "fit_error": fit_error, "convergence": convergence,
            "refined_window_um": [float(fit_lo), float(fit_hi)], "rounds": rounds,
            "warning": "A fitted far-field feature is a resonance candidate, not proof of an ideal BIC."}
