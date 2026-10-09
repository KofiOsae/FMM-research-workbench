"""Local research workbench. Run: python app.py; open http://127.0.0.1:8765/"""

from __future__ import annotations

import io
import importlib.metadata as metadata
import json
import os
import logging
import run_jobs
import usage_stats
from threading import Lock
from dataclasses import asdict
from datetime import datetime, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).parent / ".mplconfig"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from bands import BandModel, solve_bands, band_convergence, solve_band_mode, full_zone_gaps
from stack import (Layer, StackModel, solve_stack, stack_convergence, grid_convergence,
                   stack_field, field_validation, vertical_field)
from waveguide import WaveguideModel, solve_waveguide
from materials import material_catalog, import_material
from experiment import compare_or_fit
from cavity import purcell_estimate
from tmm import solve_tmm
from scattering_maps import angle_wavelength_map, kspace_map, polarization_kspace_map
from resonance import adaptive_resonance
from multi_resonance import adaptive_multi_resonance, multi_resonance_sweep
from sweep import evaluate_point, diffraction_order_sweep
from multifit import multi_parameter_fit
from slab_compare import slab_phase_match, compare_stack_layer, dispersion_comparison
from multilayer_modes import solve_multilayer_modes
from tolerance import tolerance_study
from research_report import research_report
from resonance_fields import resonance_field_report
from resonant_polarization import resonant_polarization
from leaky_modes import solve_leaky_mode
from dipole_ldos import dipole_ldos, dipole_ldos_spectrum
from polarization_winding import polarization_winding
from optimization import optimize_geometry
from bayesian import bayesian_spectrum
from constitutive import constitutive_response
from vector_modes import solve_vector_modes
from mode_coupling import (mode_port_coupling, mode_port_coupling_sweep,
                           mode_port_coupling_validation)
from coupled_branches import fit_coupled_branches
from metasurface import phase_library_and_lens
from resonator_metrics import resonator_metrics
from research_validation import (diffraction_order_map, linked_observable_sweep,
                                 settings_fingerprint, validation_report)
from finite_grating import (FiniteGratingModel, solve_finite_grating,
                            validate_finite_grating, finite_grating_spectrum,
                            finite_grating_sweep, finite_grating_tolerance,
                            optimize_finite_grating, benchmark_finite_grating)


ROOT = Path(__file__).parent
logging.basicConfig(filename=ROOT / "workbench.log", level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
MAX_BODY = 2_000_000  # imported 256×256 masks plus project metadata
PUBLIC_DEMO = os.environ.get("PUBLIC_DEMO") == "1"
_solver_gate = Lock()


def software_versions() -> dict:
    return {name: metadata.version(name) for name in ("numpy", "scipy", "grcwa", "matplotlib")}


def method_provenance(operation: str) -> dict:
    """Machine-readable method references attached to every numerical result."""
    common = [{"id": "li-1997", "title": "New formulation of the Fourier modal method for crossed surface-relief gratings",
               "doi": "10.1364/JOSAA.14.002758", "applies_to": "FMM/RCWA scattering and fields"},
              {"id": "grcwa", "title": "grcwa documentation and source",
               "url": "https://github.com/weiliangjinca/grcwa", "applies_to": "installed FMM/RCWA implementation"}]
    groups = {
        "tmm": [{"id":"yeh-1988", "title":"Optical Waves in Layered Media",
                 "applies_to":"uniform multilayer transfer matrices"}],
        "bands": [{"id":"johnson-joannopoulos-2001", "title":"Block-iterative frequency-domain methods for Maxwell equations in a planewave basis",
                   "doi":"10.1364/OE.8.000173", "applies_to":"plane-wave band eigensolver"}],
        "band_mode": [{"id":"johnson-joannopoulos-2001", "title":"Block-iterative frequency-domain methods for Maxwell equations in a planewave basis",
                       "doi":"10.1364/OE.8.000173", "applies_to":"plane-wave band eigenfield"}],
        "vector_modes": [{"id":"fallahkhair-2008", "title":"Vector finite difference modesolver for anisotropic dielectric waveguides",
                          "doi":"10.1109/JLT.2008.923643", "applies_to":"full-vector finite-difference modes"}],
        "mode_port_coupling": [{"id":"fallahkhair-2008", "title":"Vector finite difference modesolver for anisotropic dielectric waveguides",
                                "doi":"10.1109/JLT.2008.923643", "applies_to":"waveguide eigenmodes used for port overlap"},
                               {"id":"snyder-love-1983", "title":"Optical Waveguide Theory",
                                "applies_to":"mode orthogonality, excitation, and reciprocity projection"}],
        "finite_grating_coupler": [{"id":"taillaert-2006", "title":"Grating couplers for coupling between optical fibers and nanophotonic waveguides",
                                     "doi":"10.1143/JJAP.45.6071", "applies_to":"finite grating-coupler device physics"},
                                    {"id":"snyder-love-1983", "title":"Optical Waveguide Theory",
                                     "applies_to":"guided-mode normalization and reciprocity"},
                                    {"id":"oskooi-2010", "title":"MEEP: A flexible free-software package for electromagnetic simulations by the FDTD method",
                                     "doi":"10.1016/j.cpc.2009.11.008", "applies_to":"finite-domain absorbing-boundary validation context"}],
        "coupled_branches": [{"id":"hopfield-1958", "title":"Theory of the contribution of excitons to the complex dielectric constant of crystals",
                              "doi":"10.1103/PhysRev.112.1555", "applies_to":"two-oscillator branch composition"}],
        "dipole_ldos": [{"id":"novotny-hecht", "title":"Principles of Nano-Optics, planar Green tensors",
                         "applies_to":"planar electric-dipole LDOS and collection"}],
        "dipole_ldos_spectrum": [{"id":"novotny-hecht", "title":"Principles of Nano-Optics, planar Green tensors",
                                  "applies_to":"planar electric-dipole LDOS and collection"}],
        "bayesian_spectrum": [{"id":"vehtari-2021", "title":"Rank-normalization, folding, and localization: an improved R-hat",
                               "doi":"10.1214/20-BA1221", "applies_to":"diagnostic target; implementation limitations are reported"}],
        "leaky_mode": [{"id":"kristensen-2020", "title":"Modes and mode volumes of leaky optical cavities and plasmonic nanoresonators",
                        "doi":"10.1103/RevModPhys.92.031001", "applies_to":"QNM interpretation limits"}]
    }
    fmm_operations = {"solve","spectrum","field","vertical_field","angle_wavelength","kspace",
                      "polarization_kspace","resonance","multi_resonance","multi_resonance_sweep","metasurface_phase_library",
                      "sweep_point","diffraction_sweep","linked_observable_sweep","validation_report","diffraction_order_map",
                      "slab_compare","slab_dispersion","tolerance","resonance_fields",
                      "resonant_polarization","polarization_winding","optimize_geometry","measurement","multi_fit"}
    references = (common if operation in fmm_operations else []) + groups.get(operation, [])
    return {"operation": operation, "references": references,
            "note": "Citations identify numerical formulations or interpretation standards; they do not certify convergence of this result."}


def parse_stack(data: dict) -> StackModel:
    values = {**data}
    values["layers"] = tuple(Layer(**item) for item in values.get("layers", [{}]))
    return StackModel(**values)


def _spectrum_diagnostics(rows: list[dict]) -> dict:
    failed = [row for row in rows if row.get("status") not in ("converged", "unconverged")]
    grouped: dict[str, list[float]] = {}
    for row in failed:
        grouped.setdefault(str(row.get("status") or "Unknown solver failure"), []).append(
            float(row["wavelength_um"]))

    def suggestion(message: str) -> str:
        lower = message.lower()
        if "grazing" in lower or "cutoff" in lower:
            return "Move the wavelength range or incidence angle slightly so no sampled point lies exactly on the diffraction-order cutoff, then rerun."
        if (("outside" in lower and ("range" in lower or "wavelength" in lower))
                or "data cover" in lower or "data coverage" in lower):
            return "Restrict the wavelength scan to the cited material-data range or select/import optical constants covering the requested wavelengths."
        if "order budget" in lower or "fourier" in lower and "3" in lower:
            return "Use a Fourier budget from 3 to 201. After the run succeeds, increase it within that range and inspect convergence."
        if "grid size" in lower or "geometry grid" in lower:
            return "Use a geometry grid from 16 to 512 and keep at least 10–20 cells across the smallest feature."
        if "mask" in lower or "base64" in lower:
            return "Reopen the custom-mask editor, confirm nonzero width and height, reload or redraw the mask, and save the project again."
        if "singular" in lower or "eigen" in lower or "converge" in lower:
            return "Open Single wavelength at the first failed point, shift wavelength slightly, and compare higher Fourier and geometry-grid settings."
        return "Open Single wavelength at the first failed wavelength to expose the full solver error, then correct the cited material, geometry, angle, or numerical setting."

    causes = [{"message": message, "count": len(wavelengths),
               "first_wavelength_um": min(wavelengths),
               "last_wavelength_um": max(wavelengths),
               "example_wavelengths_um": wavelengths[:5],
               "suggestion": suggestion(message)}
              for message, wavelengths in grouped.items()]
    causes.sort(key=lambda item: (-item["count"], item["first_wavelength_um"]))
    return {"requested_points": len(rows),
            "computed_points": sum(row.get("status") in ("converged", "unconverged") for row in rows),
            "converged_points": sum(row.get("status") == "converged" for row in rows),
            "unconverged_points": sum(row.get("status") == "unconverged" for row in rows),
            "failed_points": len(failed), "causes": causes}


def spectrum(model: StackModel, start: float, stop: float, points: int) -> dict:
    if not np.isfinite([start, stop]).all() or not 2 <= points <= 2001 or not start < stop:
        raise ValueError("Use 2–2001 wavelengths with increasing finite limits")
    wavelengths = np.linspace(start, stop, points)
    rows = []
    uniform_stack = all(layer.kind == "uniform" for layer in model.layers)
    for point_index, wavelength in enumerate(wavelengths):
        run_jobs.progress(point_index, points, 'Wavelength spectrum')
        current = StackModel(**{**asdict(model), "wavelength_um": float(wavelength),
                                "layers": model.layers})
        try:
            if uniform_stack:
                result = solve_tmm(current)
                rows.append({"wavelength_um": float(wavelength),
                             **{key: result[key] for key in ("R", "T", "A", "r_phase_deg", "t_phase_deg")},
                             "R0": result["R"], "T0": result["T"],
                             "status": "converged", "order_change": 0.0})
                continue
            convergence = stack_convergence(current)
            result = convergence["samples"][-1]
            rows.append({"wavelength_um": float(wavelength),
                         **{key: result[key] for key in ("R", "T", "A", "R0", "T0", "r_phase_deg", "t_phase_deg")},
                         "status": "converged" if convergence["converged"]
                         and convergence["physical_balance_ok"] else "unconverged",
                         "order_change": convergence["max_change"]})
        except (ValueError, np.linalg.LinAlgError) as exc:
            rows.append({"wavelength_um": float(wavelength), "status": str(exc)})
    method = ("exact uniform-stack transfer matrix" if uniform_stack else
              "FMM with a Fourier-order check at every wavelength")
    note = ("Every layer is laterally uniform, so the exact transfer-matrix path was used. "
            "R0 and T0 equal total R and T because no diffraction orders exist." if uniform_stack else
            "Each point compares multiple Fourier orders; geometry-grid convergence is a separate check.")
    return {"rows": rows, "diagnostics": _spectrum_diagnostics(rows),
            "method": method, "note": note,
            "phase_convention": "Wrapped degrees; uniform TMM uses tangential electric amplitudes at the first/last interfaces. Patterned FMM uses the specular outgoing local s/p component matching incident polarization. Null means undefined. Power convergence does not certify phase convergence."}


def make_figure(kind: str, result: dict, title: str, fmt: str) -> bytes:
    if fmt not in ("svg", "png") or kind not in ("bands", "spectrum", "field", "vertical", "waveguide"):
        raise ValueError("Unsupported figure type")
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                         "axes.spines.right": False, "savefig.dpi": 300})
    fig, ax = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
    if kind == "spectrum":
        rows = [row for row in result["rows"] if row.get("status") == "converged"]
        if not rows:
            raise ValueError("No computed spectrum points")
        x = [row["wavelength_um"] for row in rows]
        for key, color in (("R", "#2666b4"), ("T", "#159467"), ("A", "#d9772d")):
            ax.plot(x, [row[key] for row in rows], label=key, color=color, linewidth=2)
        ax.set(xlabel="Vacuum wavelength (µm)", ylabel="Incident power fraction", ylim=(-.03, 1.03))
        ax.legend(frameon=False, ncol=3)
    elif kind == "bands":
        x = np.asarray(result["distance"])
        for name, color, style in (("TE", "#2666b4", "-"), ("TM", "#d9772d", "--")):
            for j, band in enumerate(np.asarray(result[name]).T):
                ax.plot(x, band, color=color, linestyle=style, linewidth=1.3,
                        label=name if j == 0 else None)
        ticks = result["ticks"]
        ax.set_xticks(x[ticks], result["tick_labels"])
        for tick in x[ticks[1:-1]]:
            ax.axvline(tick, color="#d0d5dd", linewidth=.8)
        ax.set(xlabel="Bloch wavevector", ylabel="Normalized frequency  a/λ")
        ax.legend(frameon=False)
    elif kind == "waveguide":
        mode = result["modes"][int(result.get("selected_mode", 0))]
        z = np.asarray(mode["z_um"])
        ax.axvspan(0, result["model"]["thickness_um"], color="#e9f3fc", label="Core")
        ax.plot(z, mode["profile"], color="#1769aa", linewidth=2,
                label=f'{mode["name"]}: n_eff={mode["n_eff"]:.6f}')
        ax.set(xlabel="Transverse z (µm)", ylabel=f'Normalized {mode["profile_quantity"]}', ylim=(-.03, 1.05))
        ax.legend(frameon=False)
    elif kind == "vertical":
        component = result.get("selected_component", "E2")
        allowed = ("E2", "H2", "Ex_abs", "Ey_abs", "Ez_abs", "Hx_abs", "Hy_abs",
                   "Hz_abs", "Sx", "Sy", "Sz", "energy_proxy", "energy_brillouin",
                   "loss_density")
        if component not in allowed:
            raise ValueError("Unknown vertical-field quantity")
        values = np.asarray(result[component], dtype=float)
        zmin, zmax = min(result["z_um"]), max(result["z_um"])
        # The model and 3D editor store finite layers from physical top to
        # bottom.  Put depth zero at the top of exported cross sections.
        extent = (0, 1, zmax, zmin)
        if component in ("Sx", "Sy", "Sz"):
            scale = max(float(np.max(np.abs(values))), 1e-12)
            image = ax.imshow(values, origin="upper", extent=extent, aspect="auto",
                              cmap="RdBu_r", vmin=-scale, vmax=scale, interpolation="nearest")
        else:
            image = ax.imshow(values, origin="upper", extent=extent, aspect="auto",
                              cmap="magma", interpolation="nearest")
        fig.colorbar(image, ax=ax, label=component)
        epsilon = np.asarray(result.get("epsilon_real", []), dtype=float)
        if epsilon.shape == values.shape and epsilon.size:
            # Draw actual material transitions from the sampled cross section.
            # A contour level is placed midway between every distinct epsilon
            # value, which supports multilayers and patterned inclusions.
            distinct = np.unique(np.round(epsilon, decimals=10))
            levels = [(a+b)/2 for a, b in zip(distinct[:-1], distinct[1:])]
            if levels:
                coordinate = np.linspace(0, 1, epsilon.shape[1])
                ax.contour(coordinate, np.asarray(result["z_um"]), epsilon, levels=levels,
                           colors="#111827", linewidths=.8, linestyles="--")
                ax.set_ylim(zmax, zmin)
        ax.set(xlabel="Position along a₁" if result["plane"] == "xz" else "Position along a₂",
               ylabel="Depth from physical top (µm)")
    else:
        component = result.get("selected_component", "intensity")
        if component not in ("intensity", "Ex_real", "Ey_real", "Ez_real"):
            raise ValueError("Unknown field component")
        values = np.asarray(result[component], dtype=float)
        if component == "intensity":
            image = ax.imshow(values, origin="lower", extent=(0, 1, 0, 1),
                              cmap="magma", interpolation="nearest")
            label = "Relative |E|²"
        else:
            scale = max(float(np.max(np.abs(values))), 1e-12)
            image = ax.imshow(values, origin="lower", extent=(0, 1, 0, 1),
                              cmap="RdBu_r", vmin=-scale, vmax=scale,
                              interpolation="nearest")
            label = "Real " + component.removesuffix("_real") + " (incident phase reference)"
        fig.colorbar(image, ax=ax, label=label)
        epsilon = np.asarray(result.get("epsilon_real", []), dtype=float)
        if epsilon.shape == values.shape and epsilon.size:
            distinct = np.unique(np.round(epsilon, decimals=10))
            levels = [(a+b)/2 for a, b in zip(distinct[:-1], distinct[1:])]
            if levels:
                ax.contour(epsilon, levels=levels, origin="lower", extent=(0, 1, 0, 1),
                           colors="#111827", linewidths=.8, linestyles="--")
        ax.set(xlabel="x / period", ylabel="y / period")
        ax.set_aspect("equal")
    ax.set_title(title)
    output = io.BytesIO()
    fig.savefig(output, format=fmt, metadata={"Creator": "FMM Research Workbench"})
    plt.close(fig)
    return output.getvalue()


class Handler(BaseHTTPRequestHandler):
    def _send(self, payload: bytes, content_type: str, status: int = 200) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(payload)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(payload)

    def do_GET(self):
        if self.path in ("/", "/index.html"):
            self._send((ROOT/"index.html").read_bytes(), "text/html; charset=utf-8")
        elif self.path in ('/workbench.js', '/workbench.css', '/FIRST_STEPS.html',
                           '/METASURFACE_GUIDE.html'):
            mime = 'text/javascript' if self.path.endswith('.js') else 'text/css' if self.path.endswith('.css') else 'text/html'
            self._send((ROOT/self.path[1:]).read_bytes(), mime+'; charset=utf-8')
        elif self.path.startswith('/api/jobs/'):
            value = run_jobs.snapshot(self.path.split('/')[-1])
            self._send(json.dumps(value or {'error':'Job expired or unknown'}).encode(), 'application/json', 200 if value else 404)
        elif self.path == '/api/job-history':
            self._send(json.dumps({'jobs': run_jobs.history()}).encode(), 'application/json')
        elif self.path == '/api/usage-stats':
            self._send(json.dumps(usage_stats.summary()).encode(), 'application/json')
        elif self.path == "/health":
            self._send(json.dumps({'status':'ok', 'root':str(ROOT.resolve()),
                'version':'scientific-workflow-2026-10-09',
                'public_demo': PUBLIC_DEMO,
                'queue': run_jobs.status_summary()}).encode(), "application/json")
        else:
            self._send(b"Not found", "text/plain", 404)

    def do_POST(self):
        if self.path == '/api/usage/visit':
            try:
                length = int(self.headers.get('Content-Length', '0'))
                if not 0 < length <= 4096:
                    raise ValueError('Request size is invalid')
                data = json.loads(self.rfile.read(length))
                result = usage_stats.record_visit(str(data.get('visitor_id', '')),
                                                  bool(data.get('new_session', True)))
                return self._send(json.dumps(result).encode(), 'application/json')
            except (ValueError, TypeError) as exc:
                return self._send(json.dumps({'error':str(exc)}).encode(), 'application/json', 400)
        if self.path.startswith('/api/forget/'):
            removed = run_jobs.forget(self.path.split('/')[-1])
            return self._send(json.dumps({'removed':removed}).encode(), 'application/json', 200 if removed else 404)
        if self.path.startswith('/api/cancel/'):
            value = run_jobs.snapshot(self.path.split('/')[-1], cancel=True)
            return self._send(json.dumps(value or {'error':'Unknown job'}).encode(), 'application/json', 200 if value else 404)
        if self.headers.get('X-Workbench-Job') == '1':
            try:
                length = int(self.headers.get('Content-Length', '0'))
                if not 0 < length <= MAX_BODY:
                    raise ValueError('Request size is invalid')
                body, path = self.rfile.read(length), self.path
                json.loads(body)
                def calculate():
                    # Reuse the synchronous dispatcher without an internal HTTP request.
                    handler = object.__new__(Handler)
                    handler.path, handler.headers, handler.rfile = path, {'Content-Length':str(len(body))}, io.BytesIO(body)
                    captured = []
                    handler._send = lambda value, mime, status=200: captured.append((json.loads(value), status))
                    with _solver_gate:
                        handler._execute_post()
                    return captured[0]
                key = run_jobs.submit(calculate, path.removeprefix('/api/'))
                return self._send(json.dumps({'job_id':key}).encode(), 'application/json', 202)
            except (ValueError, TypeError) as exc:
                return self._send(json.dumps({'error':str(exc)}).encode(), 'application/json', 400)
        if not _solver_gate.acquire(blocking=False):
            return self._send(json.dumps({'error':'Another calculation is using the solver. Submit through the Workbench queue or wait for the active run to finish.'}).encode(), 'application/json', 429)
        try:
            return self._execute_post()
        finally:
            _solver_gate.release()

    def _execute_post(self):
        if not self.path.startswith("/api/"):
            return self._send(b"Not found", "text/plain", 404)
        operation = self.path.removeprefix("/api/")
        tracked = usage_stats.tracks_calculation(operation)
        usage_started = False
        usage_finished = False
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if not 0 < length <= MAX_BODY:
                raise ValueError("Request size is invalid")
            data = json.loads(self.rfile.read(length))
            if tracked:
                usage_stats.record_calculation(operation, "started")
                usage_started = True
            if operation == "solve":
                model = parse_stack(data["model"])
                convergence = stack_convergence(model)
                high_budget = convergence["samples"][-1]["requested_budget"]
                refined = StackModel(**{**asdict(model), "layers": model.layers,
                                        "order_budget": high_budget})
                payload = {"result": solve_stack(refined), "convergence": convergence,
                           "grid_convergence": grid_convergence(refined),
                           "model": asdict(model)}
            elif operation == "spectrum":
                model = parse_stack(data["model"])
                payload = spectrum(model, float(data["start"]), float(data["stop"]), int(data["points"]))
            elif operation == "tmm":
                model = parse_stack(data["model"])
                if "start" in data:
                    start, stop, points = float(data["start"]), float(data["stop"]), int(data["points"])
                    if not np.isfinite([start, stop]).all() or not 2 <= points <= 2001 or start >= stop:
                        raise ValueError("Use 2–2001 wavelengths with increasing finite limits")
                    rows = []
                    for point_index, wavelength in enumerate(np.linspace(start, stop, points)):
                        run_jobs.progress(point_index, points, "Transfer-matrix spectrum")
                        current = StackModel(**{**asdict(model), "layers": model.layers,
                                                "wavelength_um": float(wavelength)})
                        result = solve_tmm(current)
                        rows.append({"wavelength_um": float(wavelength),
                                     **{key: result[key] for key in ("R", "T", "A", "r_phase_deg", "t_phase_deg") if key in result}})
                    payload = {"rows": rows, "method": "uniform-stack scattering recursion"}
                else:
                    payload = solve_tmm(model)
            elif operation == "field":
                model = parse_stack(data["model"])
                layer_index, z_fraction = int(data.get("layer_index", 0)), float(data.get("z_fraction", .5))
                payload = stack_field(model, layer_index, z_fraction)
                payload["validation"] = field_validation(model, layer_index, z_fraction)
            elif operation == "vertical_field":
                model = parse_stack(data["model"])
                payload = vertical_field(model, str(data.get("plane", "xz")),
                    float(data.get("fixed_fraction", .5)),
                    int(data.get("lateral_points", 81)),
                    int(data.get("points_per_layer", 25)),
                    float(data.get("exterior_depth_um", .05)))
            elif operation == "angle_wavelength":
                model = parse_stack(data["model"])
                payload = angle_wavelength_map(model, float(data["wavelength_start"]),
                    float(data["wavelength_stop"]), int(data["wavelength_points"]),
                    float(data["theta_start"]), float(data["theta_stop"]),
                    int(data["theta_points"]), str(data.get("quantity", "R")))
            elif operation == "kspace":
                model = parse_stack(data["model"])
                payload = kspace_map(model, float(data["wavelength_um"]),
                    float(data["rho_max"]), int(data["points"]),
                    str(data.get("quantity", "R")))
            elif operation == "polarization_kspace":
                model = parse_stack(data["model"])
                payload = polarization_kspace_map(model, float(data["wavelength_um"]),
                    float(data["rho_max"]), int(data["points"]),
                    str(data.get("port", "reflected")), int(data.get("order_m", 0)),
                    int(data.get("order_n", 0)), float(data.get("min_power", 1e-10)))
            elif operation == "resonance":
                model = parse_stack(data["model"])
                payload = adaptive_resonance(model, float(data["start"]), float(data["stop"]),
                    str(data.get("quantity", "R")), str(data.get("extremum", "max")),
                    int(data.get("points", 31)), int(data.get("rounds", 3)))
            elif operation == "multi_resonance":
                model = parse_stack(data["model"])
                payload = adaptive_multi_resonance(model, float(data["start"]), float(data["stop"]),
                    str(data.get("quantity", "R")), str(data.get("extremum", "both")),
                    int(data.get("global_points", 101)), int(data.get("refinement_points", 31)),
                    int(data.get("rounds", 2)), int(data.get("maximum_resonances", 4)),
                    int(data.get("background_degree", 1)))
            elif operation == "multi_resonance_sweep":
                payload = multi_resonance_sweep(parse_stack(data["model"]), str(data["parameter"]),
                    float(data["start_value"]), float(data["stop_value"]), int(data["parameter_points"]),
                    float(data["wavelength_start"]), float(data["wavelength_stop"]),
                    str(data.get("quantity", "R")), str(data.get("extremum", "both")),
                    int(data.get("global_points", 61)), int(data.get("refinement_points", 21)),
                    int(data.get("rounds", 1)), int(data.get("maximum_resonances", 4)),
                    int(data.get("background_degree", 1)))
            elif operation == "sweep_point":
                model = parse_stack(data["model"])
                payload = evaluate_point(model, str(data["x_parameter"]),
                    float(data["x"]), str(data.get("quantity", "R")),
                    str(data["y_parameter"]) if data.get("y_parameter") else None,
                    float(data["y"]) if data.get("y_parameter") else None)
            elif operation == "diffraction_sweep":
                payload = diffraction_order_sweep(parse_stack(data["model"]),
                    str(data["parameter"]), float(data["start"]), float(data["stop"]),
                    int(data["points"]))
            elif operation == "linked_observable_sweep":
                payload = linked_observable_sweep(parse_stack(data["model"]),
                    str(data["driver"]), data.get("links", []), float(data["start"]),
                    float(data["stop"]), int(data["points"]), data["observable"],
                    str(data.get("goal", "max")), bool(data.get("validate_each", False)),
                    float(data.get("tolerance", 1e-3)))
            elif operation == "validation_report":
                payload = validation_report(parse_stack(data["model"]),
                    data.get("observables", []), float(data.get("tolerance", 1e-3)),
                    bool(data.get("adaptive", True)), str(data.get("profile", "research")))
            elif operation == "diffraction_order_map":
                payload = diffraction_order_map(parse_stack(data["model"]))
            elif operation == "bands":
                band_model = BandModel(**data["model"])
                if band_model.fourier_order > 7:
                    raise ValueError("Interactive band Fourier order is capped at 7")
                payload = solve_bands(band_model)
                payload["convergence"] = band_convergence(band_model, payload)
            elif operation == "band_mode":
                band_model = BandModel(**data["model"])
                if band_model.fourier_order > 7:
                    raise ValueError("Interactive band Fourier order is capped at 7")
                payload = solve_band_mode(band_model, str(data.get("polarization", "TE")),
                    float(data.get("kx", 0)), float(data.get("ky", 0)),
                    int(data.get("band", 0)), int(data.get("resolution", 101)))
            elif operation == "full_zone_gaps":
                payload = full_zone_gaps(BandModel(**data["model"]), int(data.get("points", 9)))
            elif operation == "waveguide":
                payload = solve_waveguide(WaveguideModel(**data["model"]))
            elif operation == "slab_compare":
                stack_model = parse_stack(data["stack"])
                if "core_layer_index" in data:
                    payload = compare_stack_layer(stack_model, int(data["core_layer_index"]),
                        str(data.get("polarization", "TE")), float(data["top_n"]),
                        float(data["bottom_n"]), int(data.get("maximum_order", 2)))
                else:
                    payload = slab_phase_match(stack_model,
                        WaveguideModel(**data["slab"]), int(data.get("maximum_order", 2)))
            elif operation == "slab_dispersion":
                payload = dispersion_comparison(parse_stack(data["stack"]),
                    int(data["core_layer_index"]), str(data.get("polarization", "TE")),
                    float(data["top_n"]), float(data["bottom_n"]),
                    int(data.get("maximum_order", 2)), float(data["start"]),
                    float(data["stop"]), int(data["points"]),
                    int(data.get("mode_order", 0)))
            elif operation == "multilayer_modes":
                payload = solve_multilayer_modes(parse_stack(data["model"]),
                    str(data.get("polarization", "TE")),
                    int(data.get("grid_points", 1201)), int(data.get("modes", 8)))
            elif operation == "tolerance":
                payload = tolerance_study(parse_stack(data["model"]), data["parameters"],
                    int(data["samples"]), float(data["start"]), float(data["stop"]),
                    int(data["points"]), str(data.get("quantity", "R")),
                    str(data.get("extremum", "max")), int(data.get("seed", 12345)),
                    data.get("correlation"), int(data.get("order_m", 0)),
                    int(data.get("order_n", 0)),
                    float(data["operating_wavelength_um"]) if data.get("operating_wavelength_um") is not None else None)
            elif operation == "research_report":
                payload = {"markdown": research_report(parse_stack(data["model"]),
                    software_versions(), str(data.get("title", "Photonics simulation"))[:120])}
            elif operation == "resonance_fields":
                payload = resonance_field_report(parse_stack(data["model"]),
                    float(data["center_um"]), float(data["linewidth_um"]),
                    str(data.get("plane", "xz")), int(data.get("lateral_points", 41)),
                    int(data.get("points_per_layer", 15)))
            elif operation == "resonant_polarization":
                payload = resonant_polarization(parse_stack(data["model"]),
                    float(data["center_um"]), float(data["linewidth_um"]),
                    int(data.get("m", 0)), int(data.get("n", 0)),
                    str(data.get("port", "reflected")),
                    float(data.get("sideband_span", 3)))
            elif operation == "leaky_mode":
                payload = solve_leaky_mode(parse_stack(data["model"]),
                    float(data["center_um"]), float(data.get("initial_q", 100)))
            elif operation == "materials":
                payload = material_catalog(float(data["wavelength_um"]),
                                           str(data.get("namespace", "")))
            elif operation == "material_import":
                payload = import_material(data["name"], data["source"], data["csv"],
                                          str(data.get("namespace", "")))
            elif operation == "measurement":
                payload = compare_or_fit(parse_stack(data["model"]), data["csv"],
                                         int(data["layer_index"]), bool(data.get("fit", False)),
                                         float(data.get("lower_um", .01)),
                                         float(data.get("upper_um", 2)))
            elif operation == "multi_fit":
                payload = multi_parameter_fit(parse_stack(data["model"]), data["datasets"],
                    data["parameters"], str(data.get("method", "local")),
                    float(data.get("validation_fraction", .2)),
                    int(data.get("bootstrap", 0)), int(data.get("seed", 12345)))
            elif operation == "purcell_estimate":
                payload = purcell_estimate(**{key: float(value) for key, value in data.items()})
            elif operation == "dipole_ldos":
                payload = dipole_ldos(parse_stack(data["model"]),
                    float(data["distance_um"]), str(data.get("orientation", "isotropic")),
                    int(data.get("points", 240)), float(data.get("evanescent_limit", 30)),
                    float(data["collection_na"]) if data.get("collection_na") is not None else None)
            elif operation == "dipole_ldos_spectrum":
                payload = dipole_ldos_spectrum(parse_stack(data["model"]),
                    float(data["start_um"]), float(data["stop_um"]),
                    int(data.get("wavelength_points", 15)), float(data["distance_um"]),
                    str(data.get("orientation", "isotropic")), int(data.get("points", 120)),
                    float(data.get("evanescent_limit", 30)),
                    float(data["collection_na"]) if data.get("collection_na") is not None else None)
            elif operation == "polarization_winding":
                payload = polarization_winding(parse_stack(data["model"]),
                    float(data["center_um"]), float(data["linewidth_um"]),
                    float(data["radius"]), int(data.get("loop_points", 12)),
                    float(data.get("sideband_span", 3)), float(data.get("search_span", 2)),
                    str(data.get("quantity", "R")), str(data.get("extremum", "max")),
                    int(data.get("m", 0)), int(data.get("n", 0)),
                    str(data.get("port", "reflected")))
            elif operation == "optimize_geometry":
                payload = optimize_geometry(parse_stack(data["model"]), data["variables"],
                    data["objectives"], data.get("linear_constraints"),
                    int(data.get("generations", 4)), int(data.get("population", 5)),
                    int(data.get("seed", 12345)), bool(data.get("polish", True)),
                    int(data.get("robust_samples", 1)), float(data.get("robust_weight", .25)),
                    bool(data.get("convergence_aware", False)),
                    float(data.get("convergence_tolerance", .01)),
                    observable_definitions=data.get("observables"),
                    optical_constraints=data.get("optical_constraints"),
                    constraint_tolerance=float(data.get("constraint_tolerance", 1e-6)))
            elif operation == "bayesian_spectrum":
                payload = bayesian_spectrum(parse_stack(data["model"]), data["csv"],
                    data["parameters"], float(data["noise_sigma"]), int(data.get("draws", 1000)),
                    int(data.get("burn", 300)), int(data.get("chains", 3)), int(data.get("seed", 12345)),
                    float(data.get("noise_correlation", 0)))
            elif operation == "constitutive_response":
                payload = constitutive_response(**{key: float(value) for key, value in data.items()})
            elif operation == "vector_modes":
                payload = solve_vector_modes(**{key: (str(value) if key in ("boundary","core_shape") else
                    int(value) if key == "modes" else float(value)) for key, value in data.items()})
            elif operation == "mode_port_coupling":
                mode_inputs = {key: (str(value) if key in ("boundary", "core_shape") else
                    int(value) if key == "modes" else float(value))
                    for key, value in data["mode_solver"].items()}
                payload = mode_port_coupling(mode_inputs, data["source"])
            elif operation == "mode_port_coupling_sweep":
                mode_inputs = {key: (str(value) if key in ("boundary", "core_shape") else
                    int(value) if key == "modes" else float(value))
                    for key, value in data["mode_solver"].items()}
                payload = mode_port_coupling_sweep(mode_inputs, data["source"],
                    str(data["parameter"]), float(data["start"]), float(data["stop"]),
                    int(data["points"]))
            elif operation == "mode_port_coupling_validation":
                mode_inputs = {key: (str(value) if key in ("boundary", "core_shape") else
                    int(value) if key == "modes" else float(value))
                    for key, value in data["mode_solver"].items()}
                payload = mode_port_coupling_validation(mode_inputs, data["source"],
                    float(data.get("neff_tolerance", 5e-4)),
                    float(data.get("overlap_tolerance", 5e-3)),
                    str(data.get("profile", "research")))
            elif operation == "finite_grating_coupler":
                payload = solve_finite_grating(FiniteGratingModel(**data["model"]))
            elif operation == "finite_grating_validation":
                payload = validate_finite_grating(FiniteGratingModel(**data["model"]))
            elif operation == "finite_grating_benchmark":
                payload = benchmark_finite_grating(FiniteGratingModel(**data["model"]))
            elif operation == "finite_grating_spectrum":
                payload = finite_grating_spectrum(FiniteGratingModel(**data["model"]),
                    float(data["start"]), float(data["stop"]), int(data["points"]))
            elif operation == "finite_grating_sweep":
                payload = finite_grating_sweep(FiniteGratingModel(**data["model"]),
                    str(data["parameter"]), float(data["start"]),
                    float(data["stop"]), int(data["points"]))
            elif operation == "finite_grating_tolerance":
                payload = finite_grating_tolerance(FiniteGratingModel(**data["model"]),
                    data.get("uncertainties", []), int(data.get("samples", 20)),
                    int(data.get("seed", 12345)), float(data.get("minimum_efficiency", .5)),
                    float(data.get("maximum_reflection", .05)),
                    float(data.get("minimum_directionality", .5)))
            elif operation == "finite_grating_optimize":
                payload = optimize_finite_grating(FiniteGratingModel(**data["model"]),
                    data.get("variables", []), int(data.get("generations", 3)),
                    int(data.get("population", 5)), int(data.get("seed", 12345)),
                    float(data.get("minimum_directionality", 0)),
                    float(data.get("maximum_reflection", 1)),
                    float(data.get("maximum_power_residual", .03)),
                    int(data.get("robust_samples", 1)), data.get("uncertainty_sigma", {}),
                    float(data.get("variability_weight", 0)))
            elif operation == "coupled_branches":
                payload = fit_coupled_branches(data["parameter"], data["branch_1_um"],
                    data["branch_2_um"],
                    None if data.get("cavity_linewidth_mev") in (None, "") else float(data["cavity_linewidth_mev"]),
                    None if data.get("matter_linewidth_mev") in (None, "") else float(data["matter_linewidth_mev"]))
            elif operation == "metasurface_phase_library":
                payload = phase_library_and_lens(parse_stack(data["model"]),
                    str(data["parameter"]), float(data["start"]), float(data["stop"]),
                    int(data["points"]), str(data.get("port", "transmitted")),
                    float(data["lens_radius_um"]), float(data["focal_length_um"]),
                    int(data.get("radial_points", 101)), float(data.get("minimum_power", 0)))
            elif operation == "resonator_metrics":
                payload = resonator_metrics(float(data["wavelength_um"]),
                    float(data["linewidth_um"]), float(data["group_index"]),
                    float(data["round_trip_length_um"]),
                    None if data.get("measured_fsr_um") in (None, "") else float(data["measured_fsr_um"]),
                    None if data.get("transmission_minimum") in (None, "") else float(data["transmission_minimum"]))
            elif operation == "figure":
                fmt = data.get("format", "svg")
                figure = make_figure(data["kind"], data["result"],
                                     str(data.get("title", "Simulation result"))[:120], fmt)
                return self._send(figure, "image/svg+xml" if fmt == "svg" else "image/png")
            else:
                if usage_started:
                    usage_stats.record_calculation(operation, "failed")
                return self._send(b"Not found", "text/plain", 404)
            payload["submission"] = {"operation": operation,
                "received_utc": datetime.now(timezone.utc).isoformat(),
                "sha256": settings_fingerprint(data), "inputs": data}
            payload["software"] = software_versions()
            payload["provenance"] = method_provenance(operation)
            encoded = json.dumps(payload, allow_nan=False).encode()
            if usage_started:
                usage_stats.record_calculation(operation, "completed")
                usage_finished = True
            self._send(encoded, "application/json")
        except (ValueError, KeyError, TypeError, np.linalg.LinAlgError) as exc:
            if usage_started and not usage_finished:
                usage_stats.record_calculation(operation, "failed")
            logging.warning("Invalid request for %s: %s", self.path, exc)
            self._send(json.dumps({"error": str(exc)}).encode(), "application/json", 400)
        except Exception:
            if usage_started and not usage_finished:
                usage_stats.record_calculation(operation, "failed")
            logging.exception("Unexpected failure for %s", self.path)
            self._send(b'{"error":"Unexpected solver error. See workbench.log and the troubleshooting guide."}',
                       "application/json", 500)


class WorkbenchHTTPServer(ThreadingHTTPServer):
    daemon_threads = True
    request_queue_size = 64


def main():
    host = os.environ.get("HOST", "127.0.0.1")
    port = int(os.environ.get("PORT", "8765"))
    if not 1 <= port <= 65535:
        raise ValueError("PORT must be between 1 and 65535")
    address = (host, port)
    print(f"Open http://{address[0]}:{address[1]}/")
    WorkbenchHTTPServer(address, Handler).serve_forever()


if __name__ == "__main__":
    main()
