"""Bounded Metropolis inference for measured optical spectra."""

from __future__ import annotations
from dataclasses import replace
from run_jobs import progress

import numpy as np
from scipy.stats import norm, rankdata

from experiment import parse_measurement
from multifit import _predict, _allowed_path
from sweep import set_parameter


def _rhat(chains: np.ndarray) -> np.ndarray:
    m, n, _ = chains.shape
    means = chains.mean(axis=1); between = n*np.var(means, axis=0, ddof=1)
    within = np.mean(np.var(chains, axis=1, ddof=1), axis=0)
    variance = (n-1)/n*within+between/n
    return np.sqrt(variance/np.maximum(within, 1e-300))


def _split_rank_rhat(chains: np.ndarray) -> np.ndarray:
    """Rank-normalized split R-hat following the diagnostic in Vehtari et al."""
    m, n, p = chains.shape
    half = n//2
    if half < 4:
        return np.full(p, np.inf)
    split = np.concatenate((chains[:, :half], chains[:, -half:]), axis=0)
    ranked = np.empty_like(split, dtype=float)
    for index in range(p):
        flat = split[:, :, index].ravel()
        probabilities = (rankdata(flat, method='average')-.375)/(len(flat)+.25)
        ranked[:, :, index] = norm.ppf(np.clip(probabilities, 1e-12, 1-1e-12)).reshape(split.shape[:2])
    folded = np.abs(ranked-np.median(ranked, axis=(0, 1)))
    return np.maximum(_rhat(ranked), _rhat(folded))


def _ess(series: np.ndarray) -> float:
    x = np.asarray(series)-np.mean(series); n=len(x); variance=np.dot(x,x)/n
    if variance <= 0: return float(n)
    total=0.0
    for lag in range(1, min(n//2, 500)):
        rho=np.dot(x[:-lag],x[lag:])/(n-lag)/variance
        if rho <= 0: break
        total += rho
    return float(n/(1+2*total))


def _ess_chains(chains: np.ndarray) -> float:
    """Conservative multi-chain ESS from within-chain autocorrelation."""
    values = np.asarray(chains, dtype=float)
    return float(min(values.size, sum(_ess(row) for row in values)))


def _tail_ess(chains: np.ndarray) -> float:
    flat = np.asarray(chains).ravel()
    low, high = np.quantile(flat, [.05, .95])
    return float(min(_ess_chains((chains <= low).astype(float)),
                     _ess_chains((chains >= high).astype(float))))


def bayesian_spectrum(model, csv_text: str, parameters: list[dict], noise_sigma: float,
                      draws: int = 1000, burn: int = 300, chains: int = 3,
                      seed: int = 12345, noise_correlation: float = 0.) -> dict:
    model.validate(); wavelengths, channels, measured = parse_measurement(csv_text)
    if not 1 <= len(parameters) <= 4: raise ValueError("Use one to four inferred parameters")
    if not 200 <= draws <= 10000 or not 50 <= burn < draws or not 2 <= chains <= 4:
        raise ValueError("Use 200–10000 draws, burn-in from 50 to draws-1, and 2–4 chains")
    if not np.isfinite(noise_sigma) or not 1e-6 <= noise_sigma <= .5:
        raise ValueError("Measurement noise sigma must be from 1e-6 to 0.5 power fraction")
    if not np.isfinite(noise_correlation) or not -.95 < noise_correlation < .95:
        raise ValueError("AR(1) noise correlation must be strictly between -0.95 and 0.95")
    paths=[]; lower=[]; upper=[]; prior_mean=[]; prior_sigma=[]
    for item in parameters:
        path=str(item["path"]); lo=float(item["lower"]); hi=float(item["upper"])
        if not _allowed_path(path) or path.startswith("dataset.") or path in paths or not lo < hi:
            raise ValueError("Use unique physical parameters with increasing bounds")
        set_parameter(model,path,(lo+hi)/2); paths.append(path); lower.append(lo); upper.append(hi)
        pm=item.get("prior_mean"); ps=item.get("prior_sigma")
        if (pm is None) != (ps is None) or ps is not None and float(ps) <= 0:
            raise ValueError("A Gaussian prior requires both prior mean and positive prior sigma")
        prior_mean.append(None if pm is None else float(pm)); prior_sigma.append(None if ps is None else float(ps))
    lower=np.asarray(lower); upper=np.asarray(upper); scale=.06*(upper-lower)
    uniform=all(layer.kind=="uniform" for layer in model.layers)
    if not uniform and draws*chains*len(wavelengths)>120000:
        raise ValueError("Patterned Bayesian request is too large; reduce draws, chains, or wavelengths")
    def log_probability(x):
        if np.any(x<lower) or np.any(x>upper): return -np.inf
        current=model
        for path,value in zip(paths,x): current=set_parameter(current,path,float(value))
        prediction=_predict(current,wavelengths,channels)
        residual=(prediction-measured)/noise_sigma
        if abs(noise_correlation) > 1e-12:
            innovations = residual.copy()
            innovations[1:] = (residual[1:]-noise_correlation*residual[:-1])/np.sqrt(1-noise_correlation**2)
            value=-.5*float(np.sum(innovations*innovations))-residual.size*np.log(noise_sigma)
        else:
            value=-.5*float(np.sum(residual*residual))-residual.size*np.log(noise_sigma)
        for xi,mean,sigma in zip(x,prior_mean,prior_sigma):
            if mean is not None: value -= .5*((xi-mean)/sigma)**2+np.log(sigma)
        return value
    rng=np.random.default_rng(seed); saved=[]; acceptance=[]
    for chain in range(chains):
        chain_scale=scale.copy()
        x=lower+(upper-lower)*(chain+1)/(chains+1); lp=log_probability(x); rows=[]; accepted=0
        for iteration in range(draws):
            progress(chain*draws+iteration, chains*draws, "Bayesian draws")
            proposal=x+rng.normal(size=len(x))*chain_scale; pp=log_probability(proposal)
            if np.log(rng.random()) < pp-lp: x,lp,accepted=proposal,pp,accepted+1
            if iteration>=burn: rows.append(x.copy())
            if (iteration+1)%100==0 and iteration<burn:
                rate=accepted/(iteration+1); chain_scale*=np.exp(np.clip(rate-.28,-.1,.1))
        saved.append(rows); acceptance.append(accepted/draws)
    array=np.asarray(saved); flat=array.reshape(-1,len(paths)); rh=_split_rank_rhat(array)
    summary=[]
    for i,path in enumerate(paths):
        values=flat[:,i]; summary.append({"path":path,"mean":float(np.mean(values)),
            "std":float(np.std(values,ddof=1)),"median":float(np.median(values)),
            "credible_interval_95":[float(np.percentile(values,2.5)),float(np.percentile(values,97.5))],
            "rhat":float(rh[i]),"effective_sample_size":_ess_chains(array[:,:,i]),
            "tail_effective_sample_size":_tail_ess(array[:,:,i])})
    predictive_count=min(80 if uniform else 20, len(flat))
    predictive_indices=np.linspace(0,len(flat)-1,predictive_count,dtype=int)
    predictions=[]
    for sample_index, posterior_index in enumerate(predictive_indices):
        progress(chains*draws+sample_index, chains*draws+predictive_count, "Posterior predictive checks")
        current=model
        for path,value in zip(paths,flat[posterior_index]): current=set_parameter(current,path,float(value))
        predictions.append(_predict(current,wavelengths,channels))
    predictive=np.asarray(predictions)
    lower_prediction,median_prediction,upper_prediction=np.quantile(predictive,[.025,.5,.975],axis=0)
    covered=(measured>=lower_prediction)&(measured<=upper_prediction)
    predictive_rmse=float(np.sqrt(np.mean((median_prediction-measured)**2)))
    converged=bool(np.all(rh<1.05) and all(_ess_chains(array[:,:,i])>=100 and _tail_ess(array[:,:,i])>=50 for i in range(len(paths))) and all(.1<r<.6 for r in acceptance))
    return {"parameters":summary,"samples":{path:flat[:,i].tolist() for i,path in enumerate(paths)},
        "acceptance_fraction":acceptance,"converged":converged,"draws_per_chain":draws-burn,
        "burn_in":burn,"chains":chains,"noise_sigma":noise_sigma,"noise_correlation_ar1":noise_correlation,"channels":channels,
        "wavelength_count":len(wavelengths),"seed":seed,
        "likelihood":"Gaussian errors with user-supplied common sigma and optional within-channel AR(1) wavelength correlation",
        "posterior_predictive":{"sample_count":predictive_count,"wavelength_um":wavelengths.tolist(),"channels":channels,
            "lower_95":lower_prediction.tolist(),"median":median_prediction.tolist(),"upper_95":upper_prediction.tolist(),
            "observed":measured.tolist(),"coverage_fraction":float(np.mean(covered)),"median_rmse":predictive_rmse},
        "warning":"Posterior intervals remain conditional on the optical model, bounds, priors, and selected Gaussian noise model. Inspect rank-normalized split R-hat, bulk/tail effective sample size, trace mixing, and posterior predictive coverage."}
