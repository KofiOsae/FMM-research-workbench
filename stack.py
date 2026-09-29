"""Configurable finite-slab RCWA model; all dimensions in micrometres."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

import grcwa
import numpy as np

from materials import material_n, material_energy_terms, MATERIAL_INFO


@dataclass(frozen=True)
class Layer:
    kind: str = "rectangle"  # uniform, stripe/slot, rectangle/hole, disk/pillar, ellipse, ring
    thickness_um: float = 0.12
    background_material: str = "dielectric"  # constant n or named optical material
    background_n: float = 1.0
    feature_material: str = "gold"  # air, dielectric, silica, silicon, gold
    feature_n: float = 3.0
    fill_x: float = 0.14  # rectangle width / period, or outer radius / period
    fill_y: float = 0.40  # rectangle height / period
    inner_radius: float = 0.10  # ring inner radius / period
    offset_x: float = 0.0  # feature-center shift / first lattice coordinate
    offset_y: float = 0.0  # feature-center shift / second lattice coordinate

    def validate(self) -> None:
        numbers = [value for value in asdict(self).values() if isinstance(value, (float, int))]
        if not all(np.isfinite(value) for value in numbers):
            raise ValueError("Layer values must be finite")
        if self.kind not in ("uniform", "stripe", "slot", "rectangle", "rectangular_hole",
                             "disk", "pillar", "circular_hole", "ellipse", "ring"):
            raise ValueError("Unsupported layer shape")
        def known(value: str) -> bool:
            return value == "dielectric" or value in MATERIAL_INFO or value.startswith("custom:")
        if not known(self.feature_material) or not known(self.background_material):
            raise ValueError("Unknown layer material")
        if self.thickness_um < 0 or self.background_n <= 0 or self.feature_n <= 0:
            raise ValueError("Layer thickness cannot be negative and indices must be positive")
        if not 0 < self.fill_x < 1 or not 0 < self.fill_y < 1:
            raise ValueError("Feature fractions must be between 0 and 1")
        if not -.5 < self.offset_x < .5 or not -.5 < self.offset_y < .5:
            raise ValueError("Feature offsets must lie between -0.5 and 0.5 periods")
        if self.kind in ("disk", "pillar", "circular_hole", "ring") and self.fill_x >= 0.5:
            raise ValueError("Circular-feature radius must be less than half the period")
        if self.kind == "ring" and not 0 < self.inner_radius < self.fill_x:
            raise ValueError("Ring inner radius must lie inside the outer radius")


@dataclass(frozen=True)
class StackModel:
    wavelength_um: float = 0.905
    theta_deg: float = 0.0
    phi_deg: float = 0.0
    polarization: str = "s"
    incident_n: float = 1.5
    exit_n: float = 1.0
    period_x_um: float = 0.6
    period_y_um: float = 0.6
    lattice_angle_deg: float = 90.0
    order_budget: int = 37
    grid_size: int = 64
    layers: tuple[Layer, ...] = field(default_factory=lambda: (Layer(),))

    def validate(self) -> None:
        numbers = [value for key, value in asdict(self).items() if key != "layers" and isinstance(value, (float, int))]
        if not all(np.isfinite(value) for value in numbers):
            raise ValueError("Inputs must be finite")
        if self.polarization not in ("s", "p", "unpolarized"):
            raise ValueError("Polarization must be s, p, or unpolarized")
        if not 0 <= self.theta_deg < 89:
            raise ValueError("Incident angle must be from 0° to less than 89°")
        if not 30 <= self.lattice_angle_deg <= 150:
            raise ValueError("Lattice-vector angle must be from 30° to 150°")
        if min(self.wavelength_um, self.incident_n, self.exit_n, self.period_x_um, self.period_y_um) <= 0:
            raise ValueError("Wavelength, indices, and periods must be positive")
        if not 3 <= self.order_budget <= 101 or not 16 <= self.grid_size <= 256:
            raise ValueError("Order budget must be 3–101 and grid size 16–256")
        if not 1 <= len(self.layers) <= 64:
            raise ValueError("Use between one and 64 physical layers")
        for layer in self.layers:
            layer.validate()
        for layer in self.layers:
            material_n(layer.background_material, self.wavelength_um, layer.background_n)
            if layer.kind != "uniform":
                material_n(layer.feature_material, self.wavelength_um, layer.feature_n)


def _mask(layer: Layer, model: StackModel) -> np.ndarray:
    u = (np.arange(model.grid_size) + .5)/model.grid_size - .5
    x, y = np.meshgrid(u, u, indexing="ij")
    x = (x-layer.offset_x+.5) % 1-.5
    y = (y-layer.offset_y+.5) % 1-.5
    if layer.kind in ("stripe", "slot"):
        mask = np.abs(x) < layer.fill_x/2
    elif layer.kind in ("rectangle", "rectangular_hole"):
        mask = (np.abs(x) < layer.fill_x/2) & (np.abs(y) < layer.fill_y/2)
    elif layer.kind in ("disk", "pillar", "circular_hole"):
        mask = x*x+y*y < layer.fill_x**2
    elif layer.kind == "ellipse":
        mask = (x/(layer.fill_x/2))**2 + (y/(layer.fill_y/2))**2 < 1
    elif layer.kind == "ring":
        r2 = x*x+y*y
        mask = (r2 < layer.fill_x**2) & (r2 > layer.inner_radius**2)
    else:
        mask = np.zeros_like(x, dtype=bool)
    return mask


def _grid(layer: Layer, model: StackModel) -> np.ndarray:
    mask = _mask(layer, model)
    host = material_n(layer.background_material, model.wavelength_um, layer.background_n)
    eps = np.full(mask.shape, host**2, dtype=complex)
    if layer.kind != "uniform":
        feature_n = material_n(layer.feature_material, model.wavelength_um, layer.feature_n)
        eps[mask] = feature_n**2
    return eps


def prepare(model: StackModel):
    model.validate()
    grcwa.set_backend("numpy")
    angle = np.deg2rad(model.lattice_angle_deg)
    obj = grcwa.obj(model.order_budget, [model.period_x_um, 0],
                    [model.period_y_um*np.cos(angle), model.period_y_um*np.sin(angle)],
                    1/model.wavelength_um, np.deg2rad(model.theta_deg),
                    np.deg2rad(model.phi_deg), verbose=0)
    obj.Add_LayerUniform(0, model.incident_n**2)
    patterned = []
    for index, layer in enumerate(model.layers):
        if layer.kind == "uniform":
            obj.Add_LayerUniform(layer.thickness_um,
                                 material_n(layer.background_material, model.wavelength_um, layer.background_n)**2)
        else:
            obj.Add_LayerGrid(layer.thickness_um, model.grid_size, model.grid_size)
            patterned.append(index)
    obj.Add_LayerUniform(0, model.exit_n**2)
    obj.Init_Setup()
    for n in (model.incident_n, model.exit_n):
        kz2 = (2*np.pi*n/model.wavelength_um)**2 - obj.kx**2 - obj.ky**2
        if np.any(np.abs(kz2) < 1e-11):
            raise ValueError("An order is at grazing cutoff; adjust wavelength or angle slightly")
    if patterned:
        obj.GridLayer_geteps(np.concatenate([_grid(model.layers[i],model).ravel() for i in patterned]))
    obj.MakeExcitationPlanewave(float(model.polarization == "p"), 0,
                                 float(model.polarization == "s"), 0)
    return obj


def _port_field_fourier(obj, layer_number: int, forward, backward):
    """Cartesian Fourier fields for specified modal amplitudes at a port."""
    q = obj.q_list[layer_number]
    phi = obj.phi_list[layer_number]
    kp = obj.kp_list[layer_number]
    hxy = phi @ (forward+backward)
    hx, hy = hxy[:obj.nG], hxy[obj.nG:]
    transverse = kp @ ((phi @ ((forward-backward)/obj.omega/q)))
    ey, ex = -transverse[:obj.nG], transverse[obj.nG:]
    hz = (obj.kx*ey-obj.ky*ex)/obj.omega
    ez = (obj.ky*hx-obj.kx*hy)/obj.omega
    eps = obj.Uniform_ep_list[obj.id_list[layer_number][2]]
    ez = ez/eps
    return (ex, ey, ez), (hx, hy, hz)


def _complex_value(value) -> dict:
    value = complex(value)
    return {"real": float(value.real), "imag": float(value.imag),
            "magnitude": float(abs(value)), "phase_deg": float(np.angle(value, deg=True))}


def _channel_polarization(e_components, j: int, kx, ky, kz) -> dict:
    """Jones/Stokes data in the local propagating-order p/s basis."""
    k = np.asarray([complex(kx), complex(ky), complex(kz)])
    khat = np.real(k)/max(float(np.linalg.norm(np.real(k))), 1e-15)
    transverse = float(np.hypot(khat[0], khat[1]))
    s = np.asarray([-khat[1], khat[0], 0.0])/transverse if transverse > 1e-12 else np.asarray([0., 1., 0.])
    p = np.cross(s, khat)
    e = np.asarray([complex(component[j]) for component in e_components])
    ep, es = np.dot(p, e), np.dot(s, e)
    s0 = float(abs(ep)**2+abs(es)**2)
    s1 = float(abs(ep)**2-abs(es)**2)
    s2 = float(2*np.real(ep*np.conj(es)))
    s3 = float(-2*np.imag(ep*np.conj(es)))
    orientation = .5*np.degrees(np.arctan2(s2, s1)) if s0 else 0.0
    ellipticity = .5*np.degrees(np.arcsin(np.clip(s3/max(s0, 1e-30), -1, 1)))
    return {"basis": "local p,s", "Ep": _complex_value(ep), "Es": _complex_value(es),
            "S0": s0, "S1": s1, "S2": s2, "S3": s3,
            "normalized": {"S1": s1/max(s0, 1e-30), "S2": s2/max(s0, 1e-30),
                           "S3": s3/max(s0, 1e-30)},
            "orientation_deg": float(orientation), "ellipticity_deg": float(ellipticity)}


def diffraction_amplitudes(obj, r_power, t_power) -> list[dict]:
    """Complex Cartesian fields for reflected and transmitted diffraction orders."""
    aN, b0 = obj.GetAmplitudes(obj.Layer_N-1, 0)[0], obj.GetAmplitudes(0, 0)[1]
    zero0 = np.zeros_like(obj.a0)
    reflected_e, reflected_h = _port_field_fourier(obj, 0, zero0, b0)
    transmitted_e, transmitted_h = _port_field_fourier(obj, obj.Layer_N-1, aN, zero0)
    rows = []
    for j, (m, n) in enumerate(obj.G):
        propagating_r = bool(r_power[j] > 1e-12)
        propagating_t = bool(t_power[j] > 1e-12)
        if not (propagating_r or propagating_t):
            continue
        kr = np.sqrt(complex(obj.omega**2*obj.Uniform_ep_list[0]-obj.kx[j]**2-obj.ky[j]**2))
        kt = np.sqrt(complex(obj.omega**2*obj.Uniform_ep_list[-1]-obj.kx[j]**2-obj.ky[j]**2))
        reflected_pol = _channel_polarization(reflected_e, j, obj.kx[j], obj.ky[j], -kr)
        transmitted_pol = _channel_polarization(transmitted_e, j, obj.kx[j], obj.ky[j], kt)
        rows.append({"m": int(m), "n": int(n), "R": float(r_power[j]), "T": float(t_power[j]),
                     "reflected_angle_deg": float(np.degrees(np.arctan2(
                         np.hypot(float(np.real(obj.kx[j])), float(np.real(obj.ky[j]))), abs(float(np.real(kr)))))),
                     "transmitted_angle_deg": float(np.degrees(np.arctan2(
                         np.hypot(float(np.real(obj.kx[j])), float(np.real(obj.ky[j]))), abs(float(np.real(kt)))))),
                     "reflected": {name: _complex_value(component[j]) for name, component in
                                   zip(("Ex", "Ey", "Ez", "Hx", "Hy", "Hz"), reflected_e+reflected_h)},
                     "transmitted": {name: _complex_value(component[j]) for name, component in
                                     zip(("Ex", "Ey", "Ez", "Hx", "Hy", "Hz"), transmitted_e+transmitted_h)},
                     "reflected_polarization": reflected_pol,
                     "transmitted_polarization": transmitted_pol})
    return rows


def solve_stack(model: StackModel) -> dict:
    if model.polarization == "unpolarized":
        parts = [solve_stack(StackModel(**{**asdict(model), "layers": model.layers,
                                          "polarization": pol})) for pol in ("s", "p")]
        orders = {(item["m"], item["n"]) for part in parts for item in part["orders"]}
        lookup = [{(item["m"], item["n"]): item for item in part["orders"]} for part in parts]
        return {**{key: (parts[0][key]+parts[1][key])/2 for key in ("R", "T", "A", "R0", "T0")},
                "actual_orders": parts[0]["actual_orders"],
                "orders": [{"m": m, "n": n,
                            "R": sum(d.get((m,n), {}).get("R",0) for d in lookup)/2,
                            "T": sum(d.get((m,n), {}).get("T",0) for d in lookup)/2}
                           for m,n in sorted(orders)]}
    obj = prepare(model)
    r, t = obj.RT_Solve(normalize=1, byorder=1)
    r, t = np.asarray(r, dtype=float), np.asarray(t, dtype=float)
    R, T = float(r.sum()), float(t.sum())
    zero = np.flatnonzero(np.all(obj.G == 0, axis=1))
    i0 = int(zero[0])
    absorption = 1-R-T
    if abs(absorption) < 1e-12:
        absorption = 0.0
    amplitudes = diffraction_amplitudes(obj, r, t)
    specular = next((row for row in amplitudes if row['m']==0 and row['n']==0), None)
    component = 'Es' if model.polarization == 's' else 'Ep'
    phases = {}
    for key, port in (('r_phase_deg','reflected'), ('t_phase_deg','transmitted')):
        value = specular[port+'_polarization'][component] if specular else None
        phases[key] = value['phase_deg'] if value and value['magnitude']>1e-10 else None
    return {
        **phases,
        "R": R, "T": T, "A": absorption,
        "R0": float(r[i0]), "T0": float(t[i0]),
        "actual_orders": int(obj.nG),
        "orders": [{"m": int(m), "n": int(n), "R": float(rv), "T": float(tv)}
                   for (m,n),rv,tv in zip(obj.G,r,t) if rv > 1e-8 or tv > 1e-8],
        "order_amplitudes": amplitudes,
    }


def stack_field(model: StackModel, layer_index: int = 0, z_fraction: float = .5) -> dict:
    if model.polarization == "unpolarized":
        raise ValueError("Field maps require s or p polarization; unpolarized light has no single coherent electric field")
    if not 0 <= layer_index < len(model.layers) or not 0 <= z_fraction <= 1:
        raise ValueError("Invalid layer index or depth fraction")
    if model.layers[layer_index].kind == "uniform":
        raise ValueError("Choose a patterned layer for a horizontal field map")
    obj = prepare(model)
    e, h = obj.Solve_FieldOnGrid(layer_index+1, z_fraction*model.layers[layer_index].thickness_um)
    intensity = sum(np.abs(component)**2 for component in e)
    p = .5*np.real(np.cross(np.moveaxis(e, 0, -1), np.conj(np.moveaxis(h, 0, -1))))
    return {"intensity": intensity.T.real.tolist(),
            **{name+"_phase_deg": np.where(np.abs(component.T)>1e-10, np.angle(component.T, deg=True), None).tolist() for name, component in zip(("Ex","Ey","Ez","Hx","Hy","Hz"), (*e,*h))},
            "Ex_real": e[0].T.real.tolist(), "Ey_real": e[1].T.real.tolist(),
            "Ez_real": e[2].T.real.tolist(),
            "Hx_real": h[0].T.real.tolist(), "Hy_real": h[1].T.real.tolist(),
            "Hz_real": h[2].T.real.tolist(),
            "Sx": p[:, :, 0].T.tolist(), "Sy": p[:, :, 1].T.tolist(),
            "Sz": p[:, :, 2].T.tolist(),
            "epsilon_real": _grid(model.layers[layer_index], model).T.real.tolist(),
            "layer_index": layer_index, "z_fraction": z_fraction}


def _fourier_line(obj, layer_number: int, z_offset: float, axis: str,
                  fixed_fraction: float, samples: int):
    """Reconstruct complex fields along one primitive-cell coordinate."""
    fe, fh = obj.Solve_FieldFourier(layer_number, z_offset)
    coordinate = np.linspace(0, 1, samples, endpoint=True)
    if axis == "x":
        x = coordinate*obj.L1[0] + fixed_fraction*obj.L2[0]
        y = coordinate*obj.L1[1] + fixed_fraction*obj.L2[1]
    else:
        x = fixed_fraction*obj.L1[0] + coordinate*obj.L2[0]
        y = fixed_fraction*obj.L1[1] + coordinate*obj.L2[1]
    phase = np.exp(1j*(np.outer(obj.kx, x) + np.outer(obj.ky, y)))
    e = np.asarray([np.asarray(component) @ phase for component in fe])
    h = np.asarray([np.asarray(component) @ phase for component in fh])
    return coordinate, e, h


def _mean_normal_flux(obj, layer_number: int, z_offset: float) -> float:
    """Unit-cell average of 0.5 Re(E x H*)_z using Fourier orthogonality."""
    e, h = obj.Solve_FieldFourier(layer_number, z_offset)
    return float(.5*np.real(np.sum(e[0]*np.conj(h[1])-e[1]*np.conj(h[0]))))


def _mean_loss_density(obj, layer_number: int, z_offset: float,
                       epsilon_grid: np.ndarray, phase: np.ndarray) -> float:
    """Unit-cell average of normalized omega*Im(epsilon)*|E|²."""
    fe, _ = obj.Solve_FieldFourier(layer_number, z_offset)
    e = np.asarray([np.asarray(component) @ phase for component in fe])
    e2 = np.sum(np.abs(e)**2, axis=0).reshape(epsilon_grid.shape)
    return float(np.real(obj.normalization*obj.omega*np.mean(np.imag(epsilon_grid)*e2)))


def _section_edges(layer: Layer, plane: str, fixed_fraction: float) -> list[float]:
    """Geometric sidewall intersections in normalized primitive coordinates."""
    if layer.kind == "uniform":
        return []
    axis_x = plane == "xz"
    along_center = .5 + (layer.offset_x if axis_x else layer.offset_y)
    orth_center = .5 + (layer.offset_y if axis_x else layer.offset_x)
    distance = (fixed_fraction-orth_center+.5) % 1-.5
    half = None
    if layer.kind in ("stripe", "slot"):
        if axis_x:
            half = layer.fill_x/2
        elif abs(distance) >= layer.fill_x/2:
            return []
        else:
            return []  # the feature fills the complete displayed y direction
    elif layer.kind in ("rectangle", "rectangular_hole"):
        orth_half = layer.fill_y/2 if axis_x else layer.fill_x/2
        if abs(distance) >= orth_half:
            return []
        half = layer.fill_x/2 if axis_x else layer.fill_y/2
    elif layer.kind in ("disk", "pillar", "circular_hole", "ring"):
        if abs(distance) >= layer.fill_x:
            return []
        half = float(np.sqrt(layer.fill_x**2-distance**2))
    elif layer.kind == "ellipse":
        along = layer.fill_x/2 if axis_x else layer.fill_y/2
        orth = layer.fill_y/2 if axis_x else layer.fill_x/2
        if abs(distance) >= orth:
            return []
        half = float(along*np.sqrt(1-(distance/orth)**2))
    if half is None:
        return []
    half_widths = [half]
    if layer.kind == "ring" and abs(distance) < layer.inner_radius:
        half_widths.append(float(np.sqrt(layer.inner_radius**2-distance**2)))
    edges = []
    for width in half_widths:
        for shift in (-1, 0, 1):
            for edge in (along_center-width+shift, along_center+width+shift):
                if 1e-10 < edge < 1-1e-10:
                    edges.append(float(edge))
    return sorted(set(round(edge, 12) for edge in edges))


def vertical_field(model: StackModel, plane: str = "xz", fixed_fraction: float = .5,
                   lateral_points: int = 81, points_per_layer: int = 25,
                   exterior_depth_um: float = .05) -> dict:
    """Fields through the complete finite stack on an a1-z or a2-z section.

    Poynting quantities use 0.5 Re(E x H*) in the native grcwa normalization.
    The energy quantity is a relative nondispersive proxy, not Brillouin energy
    for a frequency-dispersive material.
    """
    if model.polarization == "unpolarized":
        raise ValueError("Vertical coherent fields require s or p polarization")
    if plane not in ("xz", "yz") or not 0 <= fixed_fraction <= 1:
        raise ValueError("Choose xz or yz and a fixed coordinate from 0 to 1")
    if not 21 <= lateral_points <= 201 or not 5 <= points_per_layer <= 81:
        raise ValueError("Use 21–201 lateral points and 5–81 points per layer")
    if not np.isfinite(exterior_depth_um) or not 0 < exterior_depth_um <= 10:
        raise ValueError("Exterior display depth must be between 0 and 10 µm")
    obj = prepare(model)
    axis = "x" if plane == "xz" else "y"
    rows = {key: [] for key in ("E2", "H2", "Ex_abs", "Ey_abs", "Ez_abs",
                                "Hx_abs", "Hy_abs", "Hz_abs", "Sx", "Sy", "Sz",
                                "Ex_real", "Ey_real", "Ez_real",
                                "Hx_real", "Hy_real", "Hz_real",
                                "energy_proxy", "energy_brillouin", "loss_density",
                                "epsilon_real")}
    rows.update({name+"_phase_deg": [] for name in ("Ex","Ey","Ez","Hx","Hy","Hz")})
    z_values, layer_numbers = [], []
    energy_materials, energy_valid = [], True
    volume_absorption = []
    integration_grid = min(model.grid_size, 96)
    u = (np.arange(integration_grid)+.5)/integration_grid-.5
    uu, vv = np.meshgrid(u, u, indexing="ij")
    gx = uu.ravel()*obj.L1[0]+vv.ravel()*obj.L2[0]
    gy = uu.ravel()*obj.L1[1]+vv.ravel()*obj.L2[1]
    integration_phase = np.exp(1j*(np.outer(obj.kx, gx)+np.outer(obj.ky, gy)))
    z_base = 0.0
    coordinate = None

    def append_sample(layer_number, z_offset, z_global, eps_line, energy_line,
                      loss_epsilon=0.0):
        nonlocal coordinate
        coordinate, e, h = _fourier_line(obj, layer_number, float(z_offset), axis,
                                          fixed_fraction, lateral_points)
        p = .5*np.real(np.cross(e.T, np.conj(h.T))).T
        e2, h2 = np.sum(np.abs(e)**2, axis=0), np.sum(np.abs(h)**2, axis=0)
        values = {"E2": e2, "H2": h2,
                  "Ex_abs": np.abs(e[0]), "Ey_abs": np.abs(e[1]), "Ez_abs": np.abs(e[2]),
                  "Hx_abs": np.abs(h[0]), "Hy_abs": np.abs(h[1]), "Hz_abs": np.abs(h[2]),
                  "Ex_real": np.real(e[0]), "Ey_real": np.real(e[1]), "Ez_real": np.real(e[2]),
                  "Hx_real": np.real(h[0]), "Hy_real": np.real(h[1]), "Hz_real": np.real(h[2]),
                  "Sx": p[0], "Sy": p[1], "Sz": p[2],
                  "energy_proxy": .25*(np.real(eps_line)*e2+h2),
                  "energy_brillouin": .25*(energy_line*e2+h2),
                  "loss_density": np.real(float(obj.normalization)*obj.omega*loss_epsilon*e2),
                  "epsilon_real": np.real(eps_line)}
        for name, component in zip(("Ex","Ey","Ez","Hx","Hy","Hz"), (*e,*h)):
            rows[name+"_phase_deg"].append(np.where(np.abs(component)>1e-10, np.angle(component, deg=True), None).tolist())
        for key, value in values.items():
            rows[key].append(np.broadcast_to(np.asarray(value, dtype=float),
                                             (lateral_points,)).tolist())
        z_values.append(float(z_global))
        layer_numbers.append(int(layer_number))

    exterior_points = max(5, points_per_layer//2)
    top_eps = np.full(lateral_points, model.incident_n**2, dtype=float)
    for z_global in np.linspace(-exterior_depth_um, 0, exterior_points, endpoint=False):
        append_sample(0, z_global, z_global, top_eps, top_eps)

    for index, layer in enumerate(model.layers):
        offsets = np.linspace(0, layer.thickness_um, points_per_layer,
                              endpoint=index == len(model.layers)-1)
        if index:
            offsets = offsets[1:]
        eps_grid = _grid(layer, model)
        integration_model = StackModel(**{**asdict(model), "layers": model.layers,
                                          "grid_size": integration_grid})
        integration_epsilon = _grid(layer, integration_model)
        host_terms = material_energy_terms(layer.background_material, model.wavelength_um,
                                           layer.background_n)
        feature_terms = host_terms if layer.kind == "uniform" else material_energy_terms(
            layer.feature_material, model.wavelength_um, layer.feature_n)
        energy_materials.append({"layer": index+1, "background": host_terms,
                                 "feature": None if layer.kind == "uniform" else feature_terms})
        energy_valid = energy_valid and host_terms["weak_loss_valid"] and feature_terms["weak_loss_valid"]
        mask = _mask(layer, model)
        energy_grid = np.full(mask.shape, host_terms["electric_energy_coefficient"], dtype=float)
        if layer.kind != "uniform":
            energy_grid[mask] = feature_terms["electric_energy_coefficient"]
        sample_index = np.minimum((np.linspace(0, 1, lateral_points)*(model.grid_size-1)).round().astype(int),
                                  model.grid_size-1)
        fixed_index = min(int(round(fixed_fraction*(model.grid_size-1))), model.grid_size-1)
        eps_line = eps_grid[sample_index, fixed_index] if axis == "x" else eps_grid[fixed_index, sample_index]
        energy_line = energy_grid[sample_index, fixed_index] if axis == "x" else energy_grid[fixed_index, sample_index]
        for offset in offsets:
            append_sample(index+1, offset, z_base+float(offset), eps_line,
                          energy_line, np.imag(eps_line))
        integration_offsets = np.linspace(0, layer.thickness_um, points_per_layer)
        loss_profile = [_mean_loss_density(obj, index+1, float(offset),
                        integration_epsilon, integration_phase) for offset in integration_offsets]
        volume_absorption.append({"layer": index+1,
            "A_volume": float(np.trapezoid(loss_profile, integration_offsets)),
            "z_points": points_per_layer, "lateral_grid": integration_grid})
        z_base += layer.thickness_um
    bottom_eps = np.full(lateral_points, model.exit_n**2, dtype=float)
    for offset in np.linspace(exterior_depth_um/exterior_points,
                              exterior_depth_um, exterior_points):
        append_sample(len(model.layers)+1, offset, z_base+offset,
                      bottom_eps, bottom_eps)
    # Flux-drop estimate at each layer's two faces, normalized to the incident flux.
    absorption = []
    # RT_Solve uses the same native flux convention and multiplies by
    # obj.normalization. The extra factor two follows grcwa's documented
    # convention (native complex-field flux carries the usual 1/2).
    flux_scale = 2*float(obj.normalization)
    for index, layer in enumerate(model.layers):
        delta = min(max(layer.thickness_um*1e-7, 1e-12), max(layer.thickness_um*.01, 1e-12))
        top_flux = _mean_normal_flux(obj, index+1, delta)
        bottom_flux = _mean_normal_flux(obj, index+1, max(layer.thickness_um-delta, 0))
        absorption.append({"layer": index+1, "top_flux": float(top_flux*flux_scale),
                           "bottom_flux": float(bottom_flux*flux_scale),
                           "A": float((top_flux-bottom_flux)*flux_scale)})
    boundaries, sidewalls, depth = [0.0], [], 0.0
    for index, layer in enumerate(model.layers):
        next_depth = depth+layer.thickness_um
        sidewalls.extend({"layer": index+1, "coordinate": edge,
                          "depth_start_um": depth, "depth_stop_um": next_depth}
                         for edge in _section_edges(layer, plane, fixed_fraction))
        boundaries.append(next_depth)
        depth = next_depth
    return {"plane": plane, "fixed_fraction": fixed_fraction,
            "coordinate": coordinate.tolist(), "z_um": z_values,
            "layer_number": layer_numbers, **rows,
            "structure_profile": {"depth_direction": "physical top to bottom",
                                  "horizontal_depths_um": boundaries,
                                  "sidewalls": sidewalls,
                                  "top_interface_um": 0.0,
                                  "bottom_interface_um": depth,
                                  "display_range_um": [-exterior_depth_um,
                                                       depth+exterior_depth_um]},
            "layer_absorption": absorption,
            "volume_absorption": volume_absorption,
            "energy_density": {"weak_loss_valid": energy_valid,
                               "materials": energy_materials},
            "notes": ["Fields are relative to the unit incident-wave convention used by grcwa.",
                      "The plotted z range includes solved fields in finite portions of both semi-infinite exterior media.",
                      "S = 0.5 Re(E x H*) in native solver units.",
                      "Weak-loss energy = 0.25{[epsilon' - lambda*d(epsilon')/dlambda]|E|^2+|H|^2}.",
                      "Loss density = normalization*omega*Im(epsilon)|E|^2 in incident-power-per-micrometre units."]}


def stack_convergence(model: StackModel, budgets=None, tolerance=.01) -> dict:
    if budgets is None:
        middle = model.order_budget
        budgets = (max(3, middle-16), middle, min(101, middle+24))
        if len(set(budgets)) != 3:
            raise ValueError("Choose an order budget that permits three distinct refinements")
    samples = []
    for budget in budgets:
        item = solve_stack(StackModel(**{**asdict(model), "order_budget": budget,
                                        "layers": model.layers}))
        samples.append({"requested_budget": budget,
                        "r_phase_deg": item.get('r_phase_deg'), "t_phase_deg": item.get('t_phase_deg'),
                        **{key: item[key] for key in ("R", "T", "A", "R0", "T0", "actual_orders")}})
    a, b = samples[-2:]
    delta = max(abs(a[key]-b[key]) for key in ("R", "T", "A"))
    distinct_orders = all(samples[i]["actual_orders"] < samples[i+1]["actual_orders"]
                          for i in range(len(samples)-1))
    physical_balance = all(np.isfinite([item[key] for key in ("R", "T", "A")]).all()
                           and all(-1e-4 <= item[key] <= 1+1e-4 for key in ("R", "T", "A"))
                           for item in samples)
    return {"samples": samples, "max_change": delta, "tolerance": tolerance,
            "distinct_orders": distinct_orders,
            "converged": bool(distinct_orders and delta <= tolerance),
            "physical_balance_ok": physical_balance}


def grid_convergence(model: StackModel, tolerance=.01) -> dict:
    """Compare the requested geometry grid with a grid twice as fine."""
    if model.grid_size > 128:
        raise ValueError("Grid convergence requires a base grid of at most 128")
    a = solve_stack(model)
    b = solve_stack(StackModel(**{**asdict(model), "layers": model.layers,
                                  "grid_size": 2*model.grid_size}))
    delta = max(abs(a[key]-b[key]) for key in ("R", "T", "A"))
    return {"grid_pair": [model.grid_size, 2*model.grid_size],
            "samples": [{"grid_size": model.grid_size, **{key: a[key] for key in ("R", "T", "A")}},
                        {"grid_size": 2*model.grid_size, **{key: b[key] for key in ("R", "T", "A")}}],
            "max_change": delta, "tolerance": tolerance,
            "converged": bool(delta <= tolerance)}


def field_validation(model: StackModel, layer_index=0, z_fraction=.5,
                     tolerance=.10) -> dict:
    """Relative RMS sensitivity of |E|² to Fourier and pixel refinement."""
    if model.grid_size > 128:
        raise ValueError("Field grid validation requires a base grid of at most 128")
    if model.order_budget >= 101:
        raise ValueError("Choose an order budget below 101 for field refinement")
    high_budget = min(101, model.order_budget+24)
    base_model = model
    high_model = StackModel(**{**asdict(model), "layers": model.layers,
                               "order_budget": high_budget})
    base = np.asarray(stack_field(base_model, layer_index, z_fraction)["intensity"])
    high = np.asarray(stack_field(high_model, layer_index, z_fraction)["intensity"])
    fine_model = StackModel(**{**asdict(high_model), "layers": model.layers,
                               "grid_size": 2*model.grid_size})
    fine = np.asarray(stack_field(fine_model, layer_index, z_fraction)["intensity"])
    fine_average = fine.reshape(model.grid_size, 2, model.grid_size, 2).mean(axis=(1,3))
    scale = max(float(np.linalg.norm(high)), 1e-12)
    order_error = float(np.linalg.norm(high-base)/scale)
    grid_error = float(np.linalg.norm(fine_average-high)/scale)
    return {"order_error_rms": order_error, "grid_error_rms": grid_error,
            "budget_pair": [model.order_budget, high_budget],
            "tolerance": tolerance,
            "converged": bool(max(order_error, grid_error) <= tolerance)}
