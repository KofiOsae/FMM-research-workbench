"""Full-vector finite-difference modes of a rectangular 2D waveguide."""

from __future__ import annotations

import os
from pathlib import Path
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).with_name(".mplconfig")))

import numpy as np
from modesolver import wgmodes


def solve_vector_modes(wavelength_um: float, core_width_um: float, core_height_um: float,
                       core_n: float, substrate_n: float, cladding_n: float,
                       padding_x_um: float = 1., padding_top_um: float = 1.,
                       padding_bottom_um: float = 1., mesh_um: float = .05,
                       modes: int = 2, guess: float | None = None,
                       boundary: str = "0000", core_shape: str = "rectangle",
                       sidewall_angle_deg: float = 90.) -> dict:
    values = np.asarray([wavelength_um, core_width_um, core_height_um, core_n,
        substrate_n, cladding_n, padding_x_um, padding_top_um, padding_bottom_um,
        mesh_um], float)
    if not np.isfinite(values).all() or min(values[:9]) <= 0 or mesh_um <= 0:
        raise ValueError("Wavelength, dimensions, indices, padding, and mesh must be positive and finite")
    if core_n <= max(substrate_n, cladding_n):
        raise ValueError("This guided-mode workspace requires core n above both surrounding indices")
    if not 1 <= modes <= 8 or boundary not in ("0000", "EEEE", "MMMM"):
        raise ValueError("Use 1–8 modes and 0000, EEEE, or MMMM boundaries")
    if core_shape not in ("rectangle", "ellipse", "trapezoid"):
        raise ValueError("Core shape must be rectangle, ellipse, or trapezoid")
    if not 20 <= sidewall_angle_deg <= 90:
        raise ValueError("Sidewall angle must be from 20° to 90°")
    width = core_width_um+2*padding_x_um
    height = padding_bottom_um+core_height_um+padding_top_um
    nx, ny = int(np.ceil(width/mesh_um)), int(np.ceil(height/mesh_um))
    if min(nx, ny) < 12 or nx*ny > 120000:
        raise ValueError(f"Requested {nx}×{ny} = {nx*ny:,} cells. Use at least 12 per direction and at most 120000 total; increase mesh spacing or reduce padding.")
    dx, dy = width/nx, height/ny
    x = (np.arange(nx)+.5)*dx-width/2
    y = (np.arange(ny)+.5)*dy-padding_bottom_um
    xx, yy = np.meshgrid(x, y)
    eps = np.where(yy < 0, substrate_n**2, cladding_n**2).astype(complex)
    if core_shape == "ellipse":
        core = (xx/(core_width_um/2))**2+((yy-core_height_um/2)/(core_height_um/2))**2 <= 1
    elif core_shape == "trapezoid":
        inset = np.maximum(0, core_height_um-yy)/np.tan(np.deg2rad(sidewall_angle_deg))
        half_width = np.maximum(0, core_width_um/2-inset)
        core = (yy >= 0) & (yy <= core_height_um) & (np.abs(xx) <= half_width)
    else:
        core = (np.abs(xx) <= core_width_um/2) & (yy >= 0) & (yy <= core_height_um)
    if np.count_nonzero(np.any(core, axis=0)) < 4 or np.count_nonzero(np.any(core, axis=1)) < 4:
        raise ValueError("Resolve the core with at least four cells across both width and height")
    eps[core] = core_n**2
    target = float(guess if guess is not None else .98*core_n)
    if not max(substrate_n, cladding_n) < target <= core_n*1.05:
        raise ValueError("The effective-index guess should lie above both surroundings and near or below core n")
    neff, ex, ey, ezj, hx, hy, hzj = wgmodes(wavelength_um, target, modes, dx, dy,
        boundary, eps=eps, collocate=True)
    neff = np.atleast_1d(neff)
    ex, ey, ezj, hx, hy, hzj = [field[..., None] if field.ndim == 2 else field
                                for field in (ex, ey, ezj, hx, hy, hzj)]
    rows = []
    for index in range(len(neff)):
        fields = [ex[:,:,index], ey[:,:,index], -1j*ezj[:,:,index],
                  hx[:,:,index], hy[:,:,index], -1j*hzj[:,:,index]]
        electric = sum(abs(item)**2 for item in fields[:3]); magnetic = sum(abs(item)**2 for item in fields[3:])
        scale = np.sqrt(max(float(np.max(electric)), 1e-300)); fields = [item/scale for item in fields]
        electric = electric/scale**2; magnetic = magnetic/scale**2
        e_transverse = float(np.sum(abs(fields[0])**2+abs(fields[1])**2))
        te_fraction = float(np.sum(abs(fields[0])**2)/max(e_transverse, 1e-300))
        core_fraction = float(np.sum(electric[core])/max(np.sum(electric), 1e-300))
        row = {"mode": index+1, "n_eff": {"real": float(np.real(neff[index])),
                "imag": float(np.imag(neff[index]))}, "te_like_fraction": te_fraction,
            "core_electric_fraction": core_fraction,
            "fields": {}, "E2": electric.real.tolist(), "H2_relative": magnetic.real.tolist()}
        for name, field in zip(("Ex","Ey","Ez","Hx","Hy","Hz"), fields):
            row["fields"][name+"_real"] = np.real(field).tolist()
            row["fields"][name+"_imag"] = np.imag(field).tolist()
            row["fields"][name+"_abs"] = np.abs(field).tolist()
            row["fields"][name+"_phase_deg"] = np.where(np.abs(field)>1e-10, np.angle(field, deg=True), np.nan).tolist()
            row["fields"][name+"_phase_deg"] = [[v if np.isfinite(v) else None for v in r] for r in row["fields"][name+"_phase_deg"]]
        rows.append(row)
    rows.sort(key=lambda row: row["n_eff"]["real"], reverse=True)
    return {"x_um": x.tolist(), "y_um": y.tolist(), "epsilon_real": eps.real.tolist(),
        "core_mask": core.tolist(), "modes": rows, "mesh": {"dx_um": dx, "dy_um": dy,
            "nx": nx, "ny": ny, "boundary": boundary},
        "geometry": {"core_width_um": core_width_um, "core_height_um": core_height_um,
            "core_shape": core_shape, "sidewall_angle_deg": sidewall_angle_deg,
            "padding_x_um": padding_x_um, "padding_top_um": padding_top_um,
            "padding_bottom_um": padding_bottom_um,
            "rasterized_core_cells_x": int(np.count_nonzero(np.any(core, axis=0))),
            "rasterized_core_cells_y": int(np.count_nonzero(np.any(core, axis=1)))},
        "method": "Fallahkhair-Li-Murphy full-vector finite difference Hx/Hy eigenproblem",
        "reference": "J. Lightwave Technol. 26, 1423-1431 (2008), DOI 10.1109/JLT.2008.923643",
        "normalization": "max(|Ex|^2+|Ey|^2+|Ez|^2)=1; field phase is arbitrary",
        "checks": ["Repeat with half the mesh spacing and larger padding.",
            "Reject modes that move materially with the outer boundary or lack core localization.",
            "For lossy/PML modes, inspect complex n_eff and boundary convergence."]}
