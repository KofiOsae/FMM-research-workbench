"""Physical regression checks. Run: .venv/Scripts/python -m unittest -v."""

import math
import base64
import unittest
import tempfile
from pathlib import Path
from unittest.mock import patch

import numpy as np
from scipy.optimize import brentq

from bands import BandModel, solve_bands, solve_band_mode
from benchmark_mpb import mpb_te_first_band
from materials import gold_n, silica_n, silicon_n, bk7_n, caf2_n, material_catalog, material_n
from stack import Layer, StackModel, solve_stack, stack_convergence, stack_field, vertical_field, _grid, _mask
from waveguide import WaveguideModel, solve_waveguide
from app import make_figure, spectrum
from tmm import solve_tmm
from experiment import compare_or_fit, parse_measurement
from cavity import purcell_estimate
from scattering_maps import angle_wavelength_map, kspace_map, polarization_kspace_map
from resonance import adaptive_resonance
from multi_resonance import fit_multi_resonance_spectrum, track_resonance_branches
from sweep import run_sweep, set_parameter
from materials import material_energy_terms
from multifit import multi_parameter_fit
from slab_compare import slab_phase_match, compare_stack_layer, dispersion_comparison
from multilayer_modes import solve_multilayer_modes
from tolerance import tolerance_study
from research_report import research_report
from resonance_fields import resonance_field_report
from resonant_polarization import resonant_polarization
from leaky_modes import solve_leaky_mode
from dipole_ldos import dipole_ldos, dipole_ldos_spectrum
from bayesian import bayesian_spectrum
from constitutive import constitutive_response
from vector_modes import solve_vector_modes
from mode_coupling import mode_port_coupling
from coupled_branches import fit_coupled_branches, HC_EV_UM
from polarization_winding import _charge_from_angles
from optimization import optimize_geometry
from resonator_metrics import resonator_metrics
from observables import evaluate_observable
from research_validation import diffraction_order_map, linked_observable_sweep, settings_fingerprint
import materials


class PhysicalValidation(unittest.TestCase):
    def test_custom_mask_resamples_with_documented_orientation(self):
        source = bytearray(8 * 8)
        source[0] = 255
        source[7] = 255
        layer = Layer(kind="custom_mask", custom_mask_width=8,
                      custom_mask_height=8,
                      custom_mask_base64=base64.b64encode(source).decode())
        model = StackModel(grid_size=16, layers=(layer,))
        mask = _mask(layer, model)
        self.assertEqual(mask.shape, (16, 16))
        self.assertTrue(mask[0, 0])
        self.assertTrue(mask[-1, 0])
        self.assertFalse(mask[0, -1])

    def test_resonator_metrics_match_definitions(self):
        result = resonator_metrics(1.55, 1.55e-6, 1.5, 3000,
                                   transmission_minimum=.25)
        self.assertAlmostEqual(result["loaded_q"], 1e6)
        self.assertAlmostEqual(result["estimated_fsr_um"], 1.55**2/(1.5*3000))
        self.assertEqual(len(result["coupling_candidates"]), 2)

    def test_uniform_spectrum_uses_exact_tmm_path(self):
        model = StackModel(wavelength_um=.55, theta_deg=45, polarization="p",
            incident_n=1.7786, exit_n=1,
            layers=(Layer(kind="uniform", thickness_um=.060,
                          background_material="silver_rakic"),
                    Layer(kind="uniform", thickness_um=.035,
                          background_material="tdbc_pva_015")))
        result = spectrum(model, .52, .62, 5)
        self.assertEqual(result["method"], "exact uniform-stack transfer matrix")
        self.assertEqual(len(result["rows"]), 5)
        for row in result["rows"]:
            self.assertEqual(row["status"], "converged")
            self.assertAlmostEqual(row["R"] + row["T"] + row["A"], 1.0, places=10)
            self.assertEqual(row["R0"], row["R"])

    def test_tdbc_pva_literature_model_produces_two_coupled_branches(self):
        at_exciton = material_n("tdbc_pva_015", 1.239841984/2.10)
        self.assertGreater(at_exciton.real, 0)
        self.assertGreater(at_exciton.imag, 0)
        model = StackModel(wavelength_um=.55, theta_deg=45,
            polarization="p", incident_n=1.7786, exit_n=1,
            layers=(
                Layer(kind="uniform", thickness_um=.060,
                      background_material="silver_rakic"),
                Layer(kind="uniform", thickness_um=.035,
                      background_material="tdbc_pva_015"),
            ))
        wavelengths = np.linspace(.47, .64, 181)
        reflectance = []
        for wavelength in wavelengths:
            current = StackModel(**{**model.__dict__, "layers": model.layers,
                                    "wavelength_um": float(wavelength)})
            reflectance.append(solve_tmm(current)["R"])
        minima = [i for i in range(1, len(reflectance)-1)
                  if reflectance[i] < reflectance[i-1]
                  and reflectance[i] < reflectance[i+1]]
        centers = wavelengths[minima]
        self.assertTrue(np.any(abs(centers-.520) < .015), centers)
        self.assertTrue(np.any(abs(centers-.606) < .015), centers)

    def test_full_vector_2d_mode_mesh_refinement(self):
        coarse = solve_vector_modes(1.55, .55, .22, 3.47, 1.44, 1,
            .825, .825, .825, .055, 2, 3.3)
        fine = solve_vector_modes(1.55, .55, .22, 3.47, 1.44, 1,
            .825, .825, .825, .0275, 2, 3.3)
        self.assertLess(abs(coarse["modes"][0]["n_eff"]["real"]-
                            fine["modes"][0]["n_eff"]["real"]), .01)
        self.assertEqual(np.asarray(fine["modes"][0]["fields"]["Ez_abs"]).shape,
                         (fine["mesh"]["ny"], fine["mesh"]["nx"]))
        self.assertGreater(fine["modes"][0]["core_electric_fraction"], .4)

    def test_mode_port_overlap_recovers_identical_imported_mode(self):
        settings = dict(wavelength_um=1.55, core_width_um=.55,
            core_height_um=.22, core_n=3.47, substrate_n=1.44,
            cladding_n=1, padding_x_um=.55, padding_top_um=.55,
            padding_bottom_um=.55, mesh_um=.055, modes=1, guess=3.3,
            boundary="0000", core_shape="rectangle", sidewall_angle_deg=90)
        solved = solve_vector_modes(**settings)
        imported = {"x_um": solved["x_um"], "y_um": solved["y_um"],
                    **solved["modes"][0]["fields"]}
        result = mode_port_coupling(settings, {"type":"imported", "field":imported})
        overlap = result["modes"][0]["power_overlap"]
        self.assertAlmostEqual(overlap["forward_efficiency"], 1, places=10)
        self.assertLess(overlap["backward_efficiency"], 1e-20)

    def test_gaussian_mode_overlap_reports_directional_port_projection(self):
        settings = dict(wavelength_um=1.55, core_width_um=.55,
            core_height_um=.22, core_n=3.47, substrate_n=1.44,
            cladding_n=1, padding_x_um=.55, padding_top_um=.55,
            padding_bottom_um=.55, mesh_um=.055, modes=2, guess=3.3,
            boundary="0000", core_shape="rectangle", sidewall_angle_deg=90)
        result = mode_port_coupling(settings, {"type":"gaussian",
            "waist_x_um":.8, "waist_y_um":.8, "center_y_um":.11,
            "polarization_angle_deg":0, "medium_n":1.44})
        self.assertTrue(result["quantity_status"]["power_overlap_available"])
        self.assertGreater(result["modes"][0]["power_overlap"]["forward_efficiency"], 0)
        self.assertLess(result["modes"][0]["power_overlap"]["forward_efficiency"], 1)

    def test_coupled_branch_fit_recovers_splitting_and_composition(self):
        parameter = np.linspace(-1, 1, 11)
        matter, coupling = 2.0, .04
        cavity = matter+.2*parameter
        root = np.sqrt(((cavity-matter)/2)**2+coupling**2)
        upper = (cavity+matter)/2+root
        lower = (cavity+matter)/2-root
        result = fit_coupled_branches(parameter, HC_EV_UM/upper,
            HC_EV_UM/lower, cavity_linewidth_mev=30, matter_linewidth_mev=20)
        self.assertAlmostEqual(result["coupling_mev"], 40, places=5)
        self.assertAlmostEqual(result["minimum_splitting_mev"], 80, places=5)
        middle = len(parameter)//2
        self.assertAlmostEqual(result["upper_photonic_fraction"][middle], .5, places=5)
        self.assertGreater(result["cooperativity_4g2_over_product"], 1)

    def test_vector_mode_curved_and_sloped_core_masks(self):
        ellipse = solve_vector_modes(1.55, .70, .40, 3.47, 1.44, 1.0,
            .60, .60, .60, .08, 1, 3.3, core_shape="ellipse")
        trapezoid = solve_vector_modes(1.55, .70, .40, 3.47, 1.44, 1.0,
            .60, .60, .60, .08, 1, 3.3, core_shape="trapezoid",
            sidewall_angle_deg=70)
        self.assertEqual(ellipse["geometry"]["core_shape"], "ellipse")
        self.assertEqual(trapezoid["geometry"]["core_shape"], "trapezoid")
        self.assertGreater(ellipse["modes"][0]["n_eff"]["real"], 1.44)
        self.assertGreater(trapezoid["modes"][0]["core_electric_fraction"], .1)

    def test_advanced_constitutive_bulk_limits(self):
        isotropic = constitutive_response(1.55, 1.5, 0, 1.5, 0,
            mu_relative=2, propagation_theta_deg=37, propagation_phi_deg=21)
        for mode in isotropic["modes"]:
            self.assertAlmostEqual(mode["n_complex"]["real"], 1.5*math.sqrt(2), places=10)
            self.assertLess(mode["longitudinal_electric_fraction"], 1e-20)
        active = constitutive_response(1.0, 1.5, 0, 1.6, 0,
            temperature_K=303.15, dn_dT_per_K=1e-5,
            intensity_W_m2=1e10, n2_m2_W=2e-20, gain_per_m=100)
        self.assertAlmostEqual(active["contributions"]["thermal_delta_n"], 1e-4)
        self.assertAlmostEqual(active["contributions"]["kerr_delta_n"], 2e-10)
        self.assertLess(active["contributions"]["gain_delta_k"], 0)

    def test_bounded_multiobjective_geometry_optimization(self):
        base = StackModel(incident_n=1, exit_n=1.5, order_budget=9, grid_size=32,
            layers=(Layer(kind="uniform", thickness_um=.14, background_n=2),))
        truth = StackModel(**{**base.__dict__, "layers": base.layers})
        targets = []
        for wavelength in (.82, 1.03):
            true_model = StackModel(**{**truth.__dict__, "layers": truth.layers,
                "wavelength_um": wavelength})
            true_model = set_parameter(true_model, "layer.0.thickness_um", .23)
            targets.append({"wavelength_um": wavelength, "quantity": "R",
                            "goal": "target", "target": solve_stack(true_model)["R"], "weight": 1})
        result = optimize_geometry(base,
            [{"path": "layer.0.thickness_um", "lower": .16, "upper": .30,
              "uncertainty_sigma": .002}],
            targets, generations=3, population=5, seed=4, polish=True,
            robust_samples=3)
        self.assertLess(abs(result["best_parameters"]["layer.0.thickness_um"]-.23), .015)
        self.assertTrue(result["validated"])
        self.assertEqual(len(result["objective_results"]), 2)
        self.assertEqual(result["robustness"]["samples"], 3)
        self.assertGreaterEqual(result["robustness"]["worst_loss"],
                                result["robustness"]["mean_loss"])

    def test_closed_loop_headless_polarization_winding(self):
        phase = np.linspace(0, 2*np.pi, 16, endpoint=False)
        charge, _ = _charge_from_angles(phase)
        self.assertAlmostEqual(charge, 1, places=12)
        zero, _ = _charge_from_angles(np.zeros(16))
        self.assertAlmostEqual(zero, 0, places=12)

    def test_planar_dipole_ldos_homogeneous_and_orientation_average(self):
        matched = StackModel(wavelength_um=.6, incident_n=1, exit_n=1,
            layers=(Layer(kind="uniform", thickness_um=.1, background_n=1),))
        result = dipole_ldos(matched, .08, points=60, evanescent_limit=15)
        self.assertAlmostEqual(result["perpendicular"], 1, places=10)
        self.assertAlmostEqual(result["parallel"], 1, places=10)
        self.assertAlmostEqual(result["collection"]["selected_collection_efficiency"], .5, places=7)
        self.assertAlmostEqual(result["isotropic"],
            (result["perpendicular"]+2*result["parallel"])/3, places=12)
        self.assertTrue(result["convergence"]["converged"])

    def test_planar_ldos_spectrum_keeps_homogeneous_normalization(self):
        matched = StackModel(wavelength_um=.6, incident_n=1.4, exit_n=1.4,
                             layers=(Layer(kind="uniform", thickness_um=.1,
                                           background_n=1.4),))
        result = dipole_ldos_spectrum(matched, .55, .65, 3, .08,
                                      points=40, evanescent_limit=10,
                                      collection_na=.9)
        self.assertEqual(len(result["rows"]), 3)
        self.assertTrue(result["all_converged"])
        for row in result["rows"]:
            self.assertAlmostEqual(row["normalized_decay_rate"], 1.0, places=6)

    def test_bayesian_uniform_film_posterior(self):
        base = StackModel(incident_n=1, exit_n=1.5, order_budget=9, grid_size=32,
            layers=(Layer(kind="uniform", thickness_um=.20, background_n=2),))
        truth = set_parameter(base, "layer.0.thickness_um", .26)
        rows = []
        for wavelength in np.linspace(.82, 1.18, 9):
            current = StackModel(**{**truth.__dict__, "layers": truth.layers,
                                    "wavelength_um": float(wavelength)})
            rows.append((wavelength, solve_tmm(current)["R"]))
        csv = "wavelength_um,R\n"+"\n".join(f"{w},{r}" for w,r in rows)
        result = bayesian_spectrum(base, csv,
            [{"path": "layer.0.thickness_um", "lower": .18, "upper": .34,
              "prior_mean": .26, "prior_sigma": .06}],
            noise_sigma=.003, draws=300, burn=100, chains=2, seed=19,
            noise_correlation=.2)
        estimate = result["parameters"][0]
        self.assertLess(abs(estimate["median"]-.26), .01)
        self.assertLess(estimate["credible_interval_95"][0], .26)
        self.assertGreater(estimate["credible_interval_95"][1], .26)
        self.assertIn("tail_effective_sample_size", estimate)
        self.assertEqual(result["posterior_predictive"]["sample_count"], 80)
        self.assertGreaterEqual(result["posterior_predictive"]["coverage_fraction"], 0)

    def test_planar_dipole_ldos_resolves_gold_near_field(self):
        gold = StackModel(wavelength_um=.6, incident_n=1, exit_n=1,
            layers=(Layer(kind="uniform", thickness_um=.05,
                          background_material="gold"),))
        result = dipole_ldos(gold, .03, points=80, evanescent_limit=30)
        self.assertGreater(result["isotropic"], 1)
        self.assertGreater(result["perpendicular_evanescent_correction"], 0)
        self.assertTrue(result["convergence"]["converged"])

    def test_correlated_tolerance_draws(self):
        model = StackModel(incident_n=1, exit_n=1.5, order_budget=9, grid_size=32,
            layers=(Layer(kind="uniform", thickness_um=.2, background_n=2),))
        result = tolerance_study(model, [
            {"path": "layer.0.thickness_um", "mean": .2, "sigma": .01},
            {"path": "layer.0.background_n", "mean": 2, "sigma": .05}],
            30, .8, 1.0, 7, correlation=[[1, .8], [.8, 1]], seed=9)
        self.assertEqual(result["input_correlation"], [[1.0, .8], [.8, 1.0]])
        self.assertGreater(result["empirical_parameter_correlation"][0][1], .45)
        with self.assertRaises(ValueError):
            tolerance_study(model, [{"path": "wavelength_um", "mean": .9, "sigma": .01}],
                5, .8, 1, 7, correlation=[[.5]])

    def test_parameter_sweep_global_and_layer_parameters(self):
        model = StackModel(incident_n=1, exit_n=1.5, order_budget=9, grid_size=32,
            layers=(Layer(kind="uniform", thickness_um=.2,
                          background_material="dielectric", background_n=2),))
        changed = set_parameter(model, "layer.0.thickness_um", .3)
        self.assertAlmostEqual(changed.layers[0].thickness_um, .3)
        result = run_sweep(model, "wavelength_um", .81, .99, 4, "R",
                           "theta_deg", 1, 11, 3)
        self.assertEqual(len(result["rows"]), 12)
        self.assertTrue(all(row["value"] is None or 0 <= row["value"] <= 1
                            for row in result["rows"]))
        self.assertIn(result["best_maximum"], result["rows"])
        with self.assertRaises(ValueError):
            set_parameter(model, "layer.2.thickness_um", .3)

    def test_signed_angle_map_and_added_feature_primitives(self):
        base = StackModel(incident_n=1, exit_n=1, period_x_um=.7, period_y_um=.7,
                          order_budget=9, grid_size=32,
                          layers=(Layer(kind="circular_hole", thickness_um=.1,
                                        background_n=2, feature_material="air", fill_x=.2),))
        angular = angle_wavelength_map(base, .9, 1.0, 3, -10, 10, 3, "R")
        self.assertAlmostEqual(angular["values"][0][1], angular["values"][2][1], places=7)
        for kind in ("slot", "rectangular_hole", "ellipse", "pillar"):
            layer = Layer(kind=kind, thickness_um=.1, background_n=2,
                          feature_material="air", fill_x=.2, fill_y=.3)
            result = solve_stack(StackModel(order_budget=9, grid_size=32, layers=(layer,)))
            self.assertTrue(np.isfinite([result["R"], result["T"], result["A"]]).all())
        centered = _grid(Layer(kind="rectangle", feature_n=2, fill_x=.2, fill_y=.2), base)
        shifted = _grid(Layer(kind="rectangle", feature_n=2, fill_x=.2, fill_y=.2,
                              offset_x=.2), base)
        self.assertFalse(np.array_equal(centered, shifted))
        with self.assertRaises(ValueError):
            Layer(offset_x=.5).validate()

    def test_adaptive_gmr_resonance_search(self):
        model = StackModel(wavelength_um=.762, polarization="s", incident_n=1,
            exit_n=1.45, period_x_um=.5, period_y_um=.5, order_budget=41, grid_size=64,
            layers=(Layer(kind="stripe", thickness_um=.15, background_material="silica",
                          feature_material="dielectric", feature_n=2.3, fill_x=.5),
                    Layer(kind="uniform", thickness_um=.30, background_material="silica")))
        result = adaptive_resonance(model, .74, .78, "R", "max", 15, 2)
        self.assertGreater(result["feature"]["R"], .95)
        self.assertTrue(.755 < result["feature"]["wavelength_um"] < .77)
        self.assertTrue(result["convergence"]["converged"])

    def test_joint_multi_resonance_fit_and_branch_tracking(self):
        wavelength = np.linspace(.70, .90, 301)
        spectrum = (.18 + .12*(wavelength-.8)
            + .42/(1+(2*(wavelength-.755)/.006)**2)
            - .31/(1+(2*(wavelength-.842)/.010)**2))
        result = fit_multi_resonance_spectrum(wavelength, spectrum, 3, "both", 1)
        components = result["best_model"]["components"]
        self.assertEqual(len(components), 2)
        centers = sorted(item["center_um"] for item in components)
        self.assertAlmostEqual(centers[0], .755, places=3)
        self.assertAlmostEqual(centers[1], .842, places=3)
        self.assertEqual(result["evidence_type"],
                         "joint phenomenological far-field line-shape fit")
        tracked = track_resonance_branches([components,
            [{**components[1], "center_um": .843},
             {**components[0], "center_um": .756}]],
            field_overlaps=[[[.05, .98], [.97, .04]]])
        self.assertEqual(len(tracked["branches"]), 2)
        self.assertAlmostEqual(tracked["branches"][0][-1]["center_um"], .756)

    def test_scattering_angle_and_kspace_maps(self):
        model = StackModel(incident_n=1, exit_n=1.5, period_x_um=.47,
                           period_y_um=.47, order_budget=9,
                           layers=(Layer(kind="uniform", thickness_um=.1,
                                         background_material="air"),))
        angular = angle_wavelength_map(model, .8, 1.0, 3, 0, 20, 3, "R")
        self.assertEqual(np.asarray(angular["values"]).shape, (3, 3))
        self.assertAlmostEqual(angular["values"][0][1], .04, places=7)
        kmap = kspace_map(model, .9, .2, 5, "R")
        values = np.asarray([[np.nan if v is None else v for v in row]
                             for row in kmap["values"]])
        self.assertAlmostEqual(values[2, 2], .04, places=7)
        self.assertTrue(np.isnan(values[0, 0]))
        self.assertAlmostEqual(values[2, 0], values[0, 2], places=7)

    def test_far_field_stokes_kspace_map_and_coherence_gate(self):
        base = StackModel(incident_n=1, exit_n=1.5, order_budget=9, grid_size=32,
            layers=(Layer(kind="uniform", thickness_um=.1,
                          background_material="air"),))
        for polarization, expected_s1 in (("s", -1), ("p", 1)):
            model = StackModel(**{**base.__dict__, "polarization": polarization})
            result = polarization_kspace_map(model, .95, .2, 5, "transmitted")
            values = np.asarray([[np.nan if value is None else value for value in row]
                                 for row in result["maps"]["S1"]])
            self.assertTrue(np.allclose(values[np.isfinite(values)], expected_s1,
                                        atol=1e-10))
            self.assertGreater(result["valid_points"], 0)
            self.assertEqual(result["stokes_convention"], "S3=-2 Im(Ep Es*)")
            component = "Es" if polarization == "s" else "Ep"
            real = result["maps"][component + "_real"][2][2]
            imag = result["maps"][component + "_imag"][2][2]
            self.assertAlmostEqual(math.hypot(real, imag), 1, places=10)
        unpolarized = StackModel(**{**base.__dict__, "polarization": "unpolarized"})
        with self.assertRaisesRegex(ValueError, "coherent"):
            polarization_kspace_map(unpolarized, .95, .2, 5)

    def test_gold_film_prism_reflectance_dip(self):
        model = StackModel(wavelength_um=.6328, polarization="p", incident_n=1.515,
                           exit_n=1.333, period_x_um=.47, period_y_um=.47,
                           layers=(Layer(kind="uniform", thickness_um=.05,
                                         background_material="gold"),))
        reflectance = {}
        for angle in (65, 72.5, 80):
            current = StackModel(**{**model.__dict__, "theta_deg": angle})
            reflectance[angle] = solve_tmm(current)["R"]
        self.assertLess(reflectance[72.5], .03)
        self.assertLess(reflectance[72.5], reflectance[65])
        self.assertLess(reflectance[72.5], reflectance[80])

    def test_purcell_ideal_and_detuned_estimate(self):
        wavelength, n, q = 1.55, 3.4, 1000
        volume = (wavelength/n)**3
        aligned = purcell_estimate(wavelength, wavelength, n, q, volume, 1)
        self.assertAlmostEqual(aligned["ideal_on_resonance"], 3*q/(4*math.pi**2))
        self.assertAlmostEqual(aligned["estimated_mode_contribution"], aligned["ideal_on_resonance"])
        detuned = purcell_estimate(wavelength, wavelength*1.01, n, q, volume, .5)
        self.assertLess(detuned["estimated_mode_contribution"], aligned["ideal_on_resonance"]*.5)
        with self.assertRaises(ValueError):
            purcell_estimate(wavelength, wavelength, n, q, volume, 1.1)

    def test_extended_dispersion_against_published_reference_values(self):
        for name, expected in (("baf2", 1.4745), ("lif", 1.3921),
                               ("kf", 1.3624), ("nacl", 1.5442)):
            value = materials.material_n(name, .5876)
            self.assertAlmostEqual(value.real, expected, places=3)
            self.assertEqual(value.imag, 0)
        with self.assertRaises(ValueError):
            materials.material_n("baf2", .1)

    def test_tio2_sapphire_and_mgf2_literature_dispersion(self):
        for name, expected in (("tio2", 2.6142), ("al2o3", 1.7682),
                               ("mgf2", 1.3777)):
            value = materials.material_n(name, .5876)
            self.assertAlmostEqual(value.real, expected, places=4)
            self.assertEqual(value.imag, 0)
        with self.assertRaises(ValueError):
            materials.material_n("tio2", 1.6)
        with self.assertRaises(ValueError):
            materials.material_n("al2o3", 5.1)

    def test_added_metals_polymers_and_silicon_nitride(self):
        for key in ("silver_rakic", "aluminum_rakic", "copper_rakic"):
            value = materials.material_n(key, .6328)
            self.assertGreater(value.real, 0)
            self.assertGreater(value.imag, 0)
        self.assertAlmostEqual(materials.material_n("pmma", .6328).real, 1.489, places=12)
        self.assertAlmostEqual(materials.material_n("polycarbonate", .833).real, 1.569, places=12)
        self.assertAlmostEqual(materials.material_n("polystyrene", .5876).real, 1.592, places=12)
        self.assertAlmostEqual(materials.material_n("zeonex", 1.052).real, 1.520, places=12)
        self.assertAlmostEqual(materials.material_n("n_sf11", .5876).real, 1.78472, places=4)
        self.assertAlmostEqual(materials.material_n("n_bak4", .5876).real, 1.56883, places=5)
        self.assertGreater(materials.material_n("sin", 1.55).real, 1.9)
        with self.assertRaises(ValueError):
            materials.material_n("pmma", 1.55)

    def test_imported_complex_material_validation_and_interpolation(self):
        with tempfile.TemporaryDirectory() as temporary:
            with patch.object(materials, "USER_DATA", Path(temporary) / "custom.json"):
                record = materials.import_material("Test film", "Test ellipsometry, 300 K",
                    "wavelength_um,n,k\n0.5,2,0.1\n1.0,2.2,0.3\n")
                self.assertEqual(record["key"], "custom:test-film")
                self.assertAlmostEqual(materials.material_n(record["key"], .75), 2.1+.2j)
                with self.assertRaises(ValueError):
                    materials.material_n(record["key"], 1.1)
                with self.assertRaises(ValueError):
                    materials.import_material("bad", "source", "wavelength_um,n,k\n.5,2,-1\n1,2,0\n")

    def test_measured_uniform_film_fit_and_invalid_units(self):
        truth = StackModel(layers=(Layer(kind="uniform", thickness_um=.26,
                      background_material="dielectric", background_n=2),))
        rows = ["wavelength_um,R,T"]
        for wavelength in np.linspace(.75, 1.15, 21):
            result = solve_tmm(StackModel(**{**truth.__dict__, "wavelength_um": float(wavelength)}))
            rows.append(f"{wavelength},{result['R']},{result['T']}")
        csv_text = "\n".join(rows)
        initial = StackModel(layers=(Layer(kind="uniform", thickness_um=.22,
                        background_material="dielectric", background_n=2),))
        result = compare_or_fit(initial, csv_text, 0, True, .15, .35)
        self.assertAlmostEqual(result["thickness_um"], .26, places=5)
        self.assertLess(result["rmse"], 1e-7)
        with self.assertRaises(ValueError):
            parse_measurement("wavelength_um,R\n500,40\n600,30\n700,20\n800,10\n900,5")

    def test_multi_parameter_fit_covariance_and_validation(self):
        truth = StackModel(incident_n=1, exit_n=1.5,
            layers=(Layer(kind="uniform", thickness_um=.26,
                          background_material="dielectric", background_n=2),))
        rows = ["wavelength_um,R,T"]
        for wavelength in np.linspace(.7, 1.3, 25):
            result = solve_tmm(StackModel(**{**truth.__dict__, "wavelength_um": float(wavelength)}))
            rows.append(f"{wavelength},{result['R']},{result['T']}")
        initial = StackModel(**{**truth.__dict__, "layers": (
            Layer(kind="uniform", thickness_um=.22,
                  background_material="dielectric", background_n=1.8),)})
        result = multi_parameter_fit(initial,
            [{"csv": "\n".join(rows), "polarization": "s"}],
            [{"path": "layer.0.thickness_um", "lower": .15, "upper": .35},
             {"path": "layer.0.background_n", "lower": 1.5, "upper": 2.5}],
            validation_fraction=.2)
        self.assertAlmostEqual(result["parameters"][0]["value"], .26, places=5)
        self.assertAlmostEqual(result["parameters"][1]["value"], 2, places=5)
        self.assertLess(result["train_rmse"], 1e-7)
        self.assertLess(result["validation_rmse"], 1e-7)
        self.assertEqual(np.asarray(result["correlation"]).shape, (2, 2))

    def test_multi_fit_recovers_wavelength_and_intensity_nuisance(self):
        model = StackModel(incident_n=1, exit_n=1.5, layers=(Layer(
            kind="uniform", thickness_um=.26, background_material="dielectric",
            background_n=2),))
        wavelengths = np.linspace(.82, 1.18, 25)
        values = [solve_tmm(StackModel(**{**model.__dict__,
                  "wavelength_um": float(w+.006)}))["R"] for w in wavelengths]
        measured = .015 + .94*np.asarray(values)
        csv = "wavelength_um,R\n" + "\n".join(
            f"{w},{r}" for w, r in zip(wavelengths, measured))
        result = multi_parameter_fit(model, [{"csv": csv, "polarization": "s"}], [
            {"path": "dataset.0.wavelength_offset_um", "lower": -.02, "upper": .02},
            {"path": "dataset.0.scale", "lower": .8, "upper": 1.1},
            {"path": "dataset.0.background", "lower": -.05, "upper": .05}],
            validation_fraction=.2)
        fitted = {item["path"]: item["value"] for item in result["parameters"]}
        self.assertAlmostEqual(fitted["dataset.0.wavelength_offset_um"], .006, places=5)
        self.assertAlmostEqual(fitted["dataset.0.scale"], .94, places=4)
        self.assertAlmostEqual(fitted["dataset.0.background"], .015, places=4)

    def test_diffraction_jones_stokes_basis_for_linear_polarizations(self):
        base = dict(incident_n=1, exit_n=1.5, theta_deg=17,
                    layers=(Layer(kind="uniform", thickness_um=.2,
                                  background_material="air"),))
        for polarization, expected_s1 in (("s", -1), ("p", 1)):
            result = solve_stack(StackModel(polarization=polarization, **base))
            order = result["order_amplitudes"][0]
            for port in ("reflected_polarization", "transmitted_polarization"):
                stokes = order[port]
                self.assertGreater(stokes["S0"], 0)
                self.assertAlmostEqual(stokes["normalized"]["S1"], expected_s1,
                                       places=10)
                self.assertAlmostEqual(stokes["normalized"]["S2"], 0, places=10)
                self.assertAlmostEqual(stokes["normalized"]["S3"], 0, places=10)
                self.assertTrue(np.isfinite([stokes["orientation_deg"],
                                             stokes["ellipticity_deg"]]).all())

    def test_resonant_polarization_sideband_decomposition(self):
        model = StackModel(polarization="s", incident_n=1, exit_n=1.5,
            layers=(Layer(kind="uniform", thickness_um=.2,
                          background_material="dielectric", background_n=2),))
        result = resonant_polarization(model, 1., .01, 0, 0, "reflected", 3)
        self.assertEqual(len(result["estimates"]), 2)
        self.assertGreaterEqual(result["sideband_stability_relative_change"], 0)
        self.assertAlmostEqual(result["driven_center"]["normalized"]["S1"], -1, places=10)
        self.assertIn("not an S-matrix pole", result["limitation"])

    def test_uniform_multilayer_complex_pole(self):
        model = StackModel(theta_deg=0, polarization="s", incident_n=1,
            exit_n=1, layers=(Layer(kind="uniform", thickness_um=.5,
                                    background_material="dielectric", background_n=2),))
        result = solve_leaky_mode(model, .67, 5)
        self.assertLess(result["total"]["denominator_residual"], 1e-7)
        self.assertGreater(result["total"]["Q"], 0)
        self.assertAlmostEqual(result["total"]["Q"], result["lossless_radiative"]["Q"], places=8)
        self.assertIsNone(result["absorptive_Q_from_inverse_Q_difference"])

    def test_transfer_matrix_matches_fmm_for_uniform_stacks(self):
        for polarization in ("s", "p", "unpolarized"):
            for angle in (0, 35, 60):
                model = StackModel(polarization=polarization, theta_deg=angle,
                    layers=(Layer(kind="uniform", background_material="gold", thickness_um=.04),
                            Layer(kind="uniform", background_material="silica", thickness_um=.2)))
                independent, fmm = solve_tmm(model), solve_stack(model)
                for quantity in ("R", "T", "A"):
                    self.assertAlmostEqual(independent[quantity], fmm[quantity], places=8)
        with self.assertRaisesRegex(ValueError, "uniform"):
            solve_tmm(StackModel(layers=(Layer(kind="stripe"),)))

    def test_johnson_christy_table_and_no_extrapolation(self):
        self.assertEqual(gold_n(.892), .17+5.663j)
        with self.assertRaises(ValueError):
            gold_n(2.0)

    def test_literature_silica_and_silicon(self):
        self.assertAlmostEqual(silica_n(.5876).real, 1.45846, places=4)
        self.assertGreater(silicon_n(.9).imag, 0)
        with self.assertRaises(ValueError):
            silicon_n(2)

    def test_added_glass_and_fluoride_dispersion(self):
        self.assertAlmostEqual(bk7_n(.5876).real, 1.5168, places=4)
        self.assertAlmostEqual(caf2_n(.5876).real, 1.4338, places=4)
        with self.assertRaises(ValueError):
            bk7_n(3)
        catalog = {row["key"]: row for row in material_catalog(.905)["materials"]}
        self.assertAlmostEqual(catalog["bk7"]["n"], bk7_n(.905).real)
        self.assertGreater(catalog["gold"]["k"], 0)

    def test_planar_guided_mode_bounds_and_field(self):
        for pol in ("TE", "TM"):
            result = solve_waveguide(WaveguideModel(polarization=pol))
            self.assertEqual(len(result["modes"]), 1)
            mode = result["modes"][0]
            self.assertGreater(mode["n_eff"], result["indices"]["bottom"])
            self.assertLess(mode["n_eff"], result["indices"]["core"])
            self.assertLess(mode["boundary_residual"], 1e-10)
            self.assertGreater(mode["confinement_fraction"], 0)
            self.assertLess(mode["confinement_fraction"], 1)
            self.assertAlmostEqual(max(mode["profile"]), 1)
            if pol == "TM":
                self.assertLess(mode["confinement_fraction"], mode["profile_integral_fraction"])
        with self.assertRaises(ValueError):
            solve_waveguide(WaveguideModel(core_material="gold"))
        self.assertEqual(solve_waveguide(WaveguideModel(thickness_um=.01,
                             core_n=1.5))["modes"], [])

    def test_guided_mode_publication_figure(self):
        result = solve_waveguide(WaveguideModel())
        self.assertIn(b"<svg", make_figure("waveguide", result, "Test mode", "svg")[:1000])

    def test_symmetric_slab_against_even_mode_equations(self):
        # Independent half-slab formulation: u tan(u) = w (TE),
        # or u tan(u) = (n_core/n_clad)^2 w (TM).
        for pol in ("TE", "TM"):
            model = WaveguideModel(wavelength_um=1, thickness_um=.5,
                        core_n=2, top_material="dielectric", top_n=1.5,
                        bottom_material="dielectric", bottom_n=1.5,
                        polarization=pol)
            mode = solve_waveguide(model)["modes"][0]
            k0 = 2*math.pi
            radius = model.thickness_um/2
            weight = 1 if pol == "TE" else (2/1.5)**2
            u_max = min(k0*radius*math.sqrt(2**2-1.5**2), math.pi/2-1e-8)
            def equation(u):
                w = math.sqrt((k0*radius)**2*(2**2-1.5**2)-u*u)
                return u*math.tan(u)-weight*w
            u = brentq(equation, 1e-9, u_max)
            expected = math.sqrt(2**2-(u/(k0*radius))**2)
            self.assertAlmostEqual(mode["n_eff"], expected, places=11)
            self.assertLess(mode["group_index_step_change"], 1e-5)

    def test_slab_grating_phase_matching_at_constructed_period(self):
        slab = WaveguideModel(wavelength_um=1, thickness_um=.5,
                              core_material="dielectric", core_n=2,
                              top_material="dielectric", top_n=1.5,
                              bottom_material="dielectric", bottom_n=1.5,
                              polarization="TE")
        beta = solve_waveguide(slab)["modes"][0]["beta_rad_um"]
        period = 2*np.pi/beta
        stack = StackModel(wavelength_um=1, period_x_um=period,
                           period_y_um=period, theta_deg=0,
                           incident_n=1, exit_n=1,
                           layers=(Layer(kind="uniform", thickness_um=.5,
                                         background_material="dielectric",
                                         background_n=2),))
        result = slab_phase_match(stack, slab, maximum_order=1)
        best = result["best_by_mode"][0]
        self.assertIn((best["m"], best["n"]), ((-1, 0), (1, 0), (0, -1), (0, 1)))
        self.assertAlmostEqual(best["relative_mismatch"], 0, places=12)
        self.assertEqual(result["equation"],
                         "beta ≈ |k_parallel + m b1 + n b2|")
        derived = compare_stack_layer(stack, 0, "TE", 1.5, 1.5, 1)
        self.assertTrue(derived["homogenization"]["exact_for_uniform_layer"])
        self.assertAlmostEqual(derived["homogenization"]["effective_core_n"], 2)
        scan = dispersion_comparison(stack, 0, "TE", 1.5, 1.5, 1,
                                     .95, 1.05, 5, 0)
        self.assertEqual(len(scan["rows"]), 5)
        self.assertIsNotNone(scan["closest_phase_match"])
        self.assertIsNotNone(scan["sampled_reflectance_peak"])

    def test_multilayer_modes_against_analytic_symmetric_slab(self):
        for polarization in ("TE", "TM"):
            analytic = solve_waveguide(WaveguideModel(wavelength_um=1,
                thickness_um=.5, core_material="dielectric", core_n=2,
                top_material="dielectric", top_n=1.5,
                bottom_material="dielectric", bottom_n=1.5,
                polarization=polarization))["modes"][0]
            stack = StackModel(wavelength_um=1, incident_n=1.5, exit_n=1.5,
                layers=(Layer(kind="uniform", thickness_um=.5,
                              background_material="dielectric", background_n=2),))
            numerical = solve_multilayer_modes(stack, polarization, 1601, 3)["modes"][0]
            self.assertAlmostEqual(numerical["n_eff"], analytic["n_eff"], delta=2e-3)
            self.assertGreater(numerical["confinement_fraction"], .4)

    def test_seeded_tolerance_study_statistics(self):
        model = StackModel(order_budget=9, grid_size=32,
            layers=(Layer(kind="uniform", thickness_um=.2,
                          background_material="dielectric", background_n=2),))
        settings = [{"path": "layer.0.thickness_um", "mean": .2,
                     "sigma": .005, "lower": .18, "upper": .22}]
        first = tolerance_study(model, settings, 5, .8, 1.0, 7, seed=7)
        second = tolerance_study(model, settings, 5, .8, 1.0, 7, seed=7)
        self.assertEqual(first["rows"], second["rows"])
        self.assertEqual(len(first["rows"]), 5)
        self.assertIn("p95", first["summary"]["feature_wavelength_um"])
        self.assertEqual(len(first["sensitivities"]), 2)

    def test_tolerance_tracks_nonzero_diffraction_order(self):
        model = StackModel(wavelength_um=.85, incident_n=1, exit_n=1,
            period_x_um=1.0, period_y_um=.4, order_budget=9, grid_size=24,
            layers=(Layer(kind="stripe", thickness_um=.2,
                          background_material="air", feature_material="dielectric",
                          feature_n=2, fill_x=.45, fill_y=.5),))
        settings = [{"path": "layer.0.thickness_um", "mean": .2,
                     "sigma": .003, "lower": .19, "upper": .21}]
        result = tolerance_study(model, settings, 5, .82, .88, 7,
            quantity="R_order", seed=11, order_m=-1, order_n=0,
            operating_wavelength_um=.85)
        self.assertEqual(result["order_m"], -1)
        self.assertIn("operating_value", result["summary"])
        self.assertTrue(all(row["feature_value"] >= 0 for row in result["rows"]))
        self.assertEqual(len(result["sensitivities"]), 3)

    def test_tolerance_rejects_mean_outside_bounds(self):
        model = StackModel(order_budget=9, grid_size=24,
            layers=(Layer(kind="uniform", thickness_um=.2,
                          background_material="dielectric", background_n=2),))
        with self.assertRaisesRegex(ValueError, "mean 1.4 lies outside its bounds 0.12–0.18"):
            tolerance_study(model, [{"path":"period_x_um", "mean":1.4,
                "sigma":.005, "lower":.12, "upper":.18}], 5, .8, 1.0, 7)

    def test_tolerance_explains_all_zero_diffraction_order(self):
        model = StackModel(wavelength_um=1.4, incident_n=1, exit_n=1,
            period_x_um=.18, period_y_um=.18, order_budget=9, grid_size=24,
            layers=(Layer(kind="stripe", thickness_um=.2,
                          background_material="air", feature_material="dielectric",
                          feature_n=2, fill_x=.45, fill_y=.5),))
        with self.assertRaisesRegex(ValueError, r"Selected order \(-1,0\) has zero"):
            tolerance_study(model, [{"path":"period_x_um", "mean":.18,
                "sigma":.001, "lower":.17, "upper":.19}], 5, 1.3, 1.45, 7,
                quantity="R_order", order_m=-1, order_n=0)

    def test_research_report_contains_reproducibility_fields(self):
        model = StackModel(layers=(Layer(kind="uniform", thickness_um=.2,
                                         background_material="silica"),))
        report = research_report(model, {"numpy": "test", "grcwa": "test"}, "Study")
        for text in ("# Study", "Material provenance", "Requested FMM harmonics",
                     "Geometry-grid refinement", "Interpretation limits", "Fused silica"):
            self.assertIn(text, report)

    def test_resonance_field_report_samples_expected_wavelengths(self):
        model = StackModel(wavelength_um=.9, order_budget=9, grid_size=32,
            layers=(Layer(kind="uniform", thickness_um=.2,
                          background_material="dielectric", background_n=2),))
        report = resonance_field_report(model, .91, .01, "xz", 21, 5)
        self.assertEqual(len(report["cases"]), 4)
        self.assertAlmostEqual(report["cases"][1]["wavelength_um"], .91)
        self.assertTrue(np.isfinite(report["center_to_off_max_E2"]))
        self.assertTrue(all(abs(case["R"]+case["T"]+case["A"]-1) < 1e-8
                            for case in report["cases"]))

    def test_stripe_is_invariant_along_y(self):
        model = StackModel(grid_size=32, layers=(Layer(kind="stripe", fill_x=.4,
                            feature_material="dielectric", feature_n=2),))
        eps = _grid(model.layers[0], model)
        self.assertTrue(np.all(eps == eps[:, :1]))
        result = solve_stack(model)
        self.assertAlmostEqual(result["R"] + result["T"], 1, places=8)

    def test_oblique_hexagonal_fmm_cell_conserves_lossless_power(self):
        model = StackModel(lattice_angle_deg=60, period_x_um=.7, period_y_um=.7,
            order_budget=19, grid_size=48,
            layers=(Layer(kind="disk", thickness_um=.12,
                          background_material="air", feature_material="dielectric",
                          feature_n=2, fill_x=.2),))
        result = solve_stack(model)
        self.assertAlmostEqual(result["R"] + result["T"], 1, places=7)
        self.assertAlmostEqual(result["A"], 0, places=7)

    def test_unpolarized_is_incoherent_mean(self):
        def run(pol):
            return solve_stack(StackModel(theta_deg=30, polarization=pol,
                         layers=(Layer(kind="uniform"),)))
        s, p, u = (run(pol) for pol in ("s", "p", "unpolarized"))
        for key in ("R", "T", "A", "R0", "T0"):
            self.assertAlmostEqual(u[key], (s[key]+p[key])/2, places=10)

    def test_fresnel_interface_normal_and_oblique(self):
        for angle in (0, 30):
            ci = math.cos(math.radians(angle))
            ct = math.sqrt(1-(1.5*math.sin(math.radians(angle)))**2)
            for pol in ("s", "p"):
                if pol == "s":
                    expected = ((1.5*ci-ct)/(1.5*ci+ct))**2
                else:
                    expected = ((ci-1.5*ct)/(ci+1.5*ct))**2
                result = solve_stack(StackModel(theta_deg=angle, polarization=pol,
                               layers=(Layer(kind="uniform", background_n=1),)))
                self.assertAlmostEqual(result["R"], expected, places=10)
                self.assertAlmostEqual(result["R"]+result["T"], 1, places=10)
                if angle == 0 and pol == "s":
                    order = result["order_amplitudes"][0]
                    self.assertGreater(order["reflected"]["Ey"]["magnitude"], 0)
                    self.assertTrue(np.isfinite(order["transmitted"]["Ey"]["phase_deg"]))

    def test_uniform_bk7_film_against_fabry_perot_formula(self):
        wavelength, thickness, n1 = .5876, .2, bk7_n(.5876).real
        model = StackModel(wavelength_um=wavelength, incident_n=1, exit_n=1,
                           layers=(Layer(kind="uniform", thickness_um=thickness,
                                         background_material="bk7"),))
        result = solve_stack(model)
        r01, r12 = (1-n1)/(1+n1), (n1-1)/(n1+1)
        phase = np.exp(2j*2*np.pi*n1*thickness/wavelength)
        expected_r = abs((r01+r12*phase)/(1+r01*r12*phase))**2
        self.assertAlmostEqual(result["R"], expected_r, places=10)
        self.assertAlmostEqual(result["R"]+result["T"], 1, places=10)

    def test_reversing_stack_preserves_reciprocal_transmission(self):
        first = Layer(kind="uniform", thickness_um=.15, background_n=2)
        second = Layer(kind="uniform", thickness_um=.23, background_n=1.3)
        forward = solve_stack(StackModel(wavelength_um=.8, incident_n=1.5,
                              exit_n=1, layers=(first,second)))
        reverse = solve_stack(StackModel(wavelength_um=.8, incident_n=1,
                              exit_n=1.5, layers=(second,first)))
        self.assertAlmostEqual(forward["T"], reverse["T"], places=10)

    def test_absorbing_gold_film_against_complex_fabry_perot(self):
        wavelength, thickness = .905, .05
        n1 = gold_n(wavelength)
        phase = np.exp(1j*2*np.pi*n1*thickness/wavelength)
        r01, r12 = (1-n1)/(1+n1), (n1-1)/(n1+1)
        denominator = 1+r01*r12*phase**2
        expected_r = abs((r01+r12*phase**2)/denominator)**2
        expected_t = abs((2/(1+n1))*(2*n1/(n1+1))*phase/denominator)**2
        result = solve_stack(StackModel(wavelength_um=wavelength,
                          incident_n=1, exit_n=1,
                          layers=(Layer(kind="uniform", thickness_um=thickness,
                                        background_material="gold"),)))
        self.assertAlmostEqual(result["R"], expected_r, places=10)
        self.assertAlmostEqual(result["T"], expected_t, places=10)
        self.assertGreater(result["A"], 0)

    def test_dielectric_pattern_conserves_power_and_field_is_finite(self):
        model = StackModel(order_budget=21, grid_size=32,
                           layers=(Layer(kind="disk", background_n=1,
                                         feature_material="dielectric", feature_n=2,
                                         fill_x=.2),))
        result = solve_stack(model)
        self.assertLess(abs(result["R"]+result["T"]-1), 1e-8)
        field = stack_field(model)
        self.assertEqual(np.asarray(field["intensity"]).shape, (32,32))
        self.assertTrue(np.isfinite(field["intensity"]).all())
        self.assertTrue(np.isfinite(field["Sx"]).all())
        self.assertEqual(np.asarray(field["Hx_real"]).shape, (32,32))
        self.assertIn(b"<svg", make_figure("field", {**field,
                      "selected_component": "Ex_real"}, "Real field", "svg")[:1000])

    def test_vertical_fields_poynting_and_layer_absorption(self):
        model = StackModel(wavelength_um=.905, incident_n=1, exit_n=1,
            order_budget=9, grid_size=32,
            layers=(Layer(kind="uniform", thickness_um=.05,
                          background_material="gold"),))
        result = vertical_field(model, "xz", .5, 31, 7)
        self.assertLess(result["z_um"][0], 0)
        self.assertGreater(result["z_um"][-1], sum(layer.thickness_um for layer in model.layers))
        self.assertEqual(result["structure_profile"]["top_interface_um"], 0)
        self.assertEqual(np.asarray(result["epsilon_real"]).shape,
                         np.asarray(result["E2"]).shape)
        patterned = StackModel(wavelength_um=.905, incident_n=1, exit_n=1,
            order_budget=9, grid_size=32,
            layers=(Layer(kind="stripe", thickness_um=.05,
                          background_material="air", feature_material="gold",
                          fill_x=.45),))
        patterned_result = vertical_field(patterned, "xz", .5, 31, 7)
        self.assertGreater(np.ptp(np.asarray(patterned_result["epsilon_real"])), 0)
        self.assertEqual(patterned_result["structure_profile"]["horizontal_depths_um"],
                         [0.0, .05])
        self.assertEqual(len(patterned_result["structure_profile"]["sidewalls"]), 2)
        two_layer = StackModel(wavelength_um=.905, incident_n=1, exit_n=1.45,
            order_budget=9, grid_size=32,
            layers=(Layer(kind="uniform", thickness_um=.11, background_material="air"),
                    Layer(kind="stripe", thickness_um=.0387, background_material="air",
                          feature_material="dielectric", feature_n=2, fill_x=.45)))
        profile = vertical_field(two_layer, "xz", .5, 31, 7)["structure_profile"]
        self.assertEqual(profile["horizontal_depths_um"], [0.0, .11, .1487])
        self.assertEqual([wall["coordinate"] for wall in profile["sidewalls"]],
                         [.275, .725])
        self.assertEqual(np.asarray(result["E2"]).shape,
                         (len(result["z_um"]), 31))
        for key in ("E2", "H2", "Ex_real", "Hz_real", "Sx", "Sy", "Sz", "energy_proxy",
                    "energy_brillouin", "loss_density"):
            self.assertTrue(np.isfinite(result[key]).all())
        self.assertFalse(result["energy_density"]["weak_loss_valid"])
        resolved = vertical_field(model, "xz", .5, 31, 41)
        line_loss = np.mean(np.asarray(resolved["loss_density"]), axis=1)
        finite = np.asarray(resolved["layer_number"]) == 1
        self.assertAlmostEqual(np.trapezoid(line_loss[finite],
                                            np.asarray(resolved["z_um"])[finite]),
                               solve_stack(model)["A"], places=3)
        self.assertAlmostEqual(resolved["volume_absorption"][0]["A_volume"],
                               resolved["layer_absorption"][0]["A"], places=3)
        self.assertAlmostEqual(sum(row["A"] for row in result["layer_absorption"]),
                               solve_stack(model)["A"], places=7)
        figure = make_figure("vertical", {**result, "selected_component": "Sz"},
                             "Vertical power flow", "svg")
        self.assertIn(b"<svg", figure[:1000])
        silica = material_energy_terms("silica", .905)
        self.assertTrue(silica["weak_loss_valid"])
        self.assertGreater(silica["electric_energy_coefficient"], silica["epsilon_real"])

    def test_uniform_order_convergence(self):
        report = stack_convergence(StackModel(layers=(Layer(kind="uniform", background_n=1),)))
        self.assertTrue(report["converged"])
        self.assertTrue(report["physical_balance_ok"])

    def test_gold_result_is_rejected_by_order_gate(self):
        # This is a known difficult case; a plausible energy balance is insufficient.
        report = stack_convergence(StackModel())
        self.assertFalse(report["converged"])

    def test_uniform_band_matches_folded_plane_waves(self):
        result = solve_bands(BandModel(background_n=2, inclusion_n=2,
                              fourier_order=1, points_per_segment=3, bands=3))
        self.assertAlmostEqual(result["TE"][0][0], 0, places=10)
        self.assertAlmostEqual(result["TM"][0][1], .5, places=10)
        self.assertAlmostEqual(result["TE"][3][0], .25, places=10)
        self.assertAlmostEqual(result["TM"][3][0], .25, places=10)

    def test_uniform_band_mode_reconstruction_and_parity(self):
        model = BandModel(background_n=2, inclusion_n=2, fourier_order=1,
                          points_per_segment=3, bands=3, grid_size=64)
        mode = solve_band_mode(model, "TM", 0, 0, 0, 65)
        self.assertLess(np.std(mode["intensity"]), 1e-10)
        self.assertAlmostEqual(mode["normalized_frequency"], 0, places=10)
        self.assertGreater(mode["inversion_overlap"]["real"], .999999)
        off_symmetry = solve_band_mode(model, "TE", .13, .07, 0, 33)
        self.assertIsNone(off_symmetry["inversion_overlap"])

    def test_band_material_contrast_changes_te_tm(self):
        result = solve_bands(BandModel(fourier_order=2, points_per_segment=3, bands=4))
        self.assertGreater(np.max(np.abs(np.asarray(result["TE"])-np.asarray(result["TM"]))), .01)

    def test_rectangular_and_hexagonal_band_lattices(self):
        for lattice, labels in (("rectangular", ["Γ", "X", "M", "Γ"]),
                                ("hexagonal", ["Γ", "M", "K", "Γ"])):
            result = solve_bands(BandModel(lattice=lattice, aspect_ratio=1.4,
                fourier_order=2, points_per_segment=3, bands=3, grid_size=64))
            self.assertEqual(result["tick_labels"], labels)
            self.assertTrue(np.isfinite(result["TE"]).all())
            self.assertTrue(np.isfinite(result["TM"]).all())
            self.assertTrue(np.all(np.diff(result["distance"]) >= 0))

    def test_mpb_square_rod_reference(self):
        # MPB tutorial: epsilon=12, r/a=.2, k=(.3,.3), first TE band=0.372604.
        self.assertLess(abs(mpb_te_first_band(7)-.372604), .0031)


class ObservableTrustTests(unittest.TestCase):
    def test_custom_diffraction_formula_and_imbalance(self):
        result = {"R": .1, "T": .9, "A": 0, "R0": .1, "T0": .2,
                  "orders": [{"m": 1, "n": 0, "R": 0, "T": .35},
                             {"m": -1, "n": 0, "R": 0, "T": .35}]}
        self.assertAlmostEqual(evaluate_observable(result, "To(1,0)+To(-1,0)"), .7)
        self.assertAlmostEqual(evaluate_observable(result, "T-(To(1,0)+To(-1,0))"), .2)
        self.assertAlmostEqual(evaluate_observable(result, "stdev([To(1,0),To(-1,0)])"), 0)
        with self.assertRaises(ValueError):
            evaluate_observable(result, "__import__('os').system('echo no')")

    def test_linked_sweep_records_relationship_and_fingerprint(self):
        model = StackModel(layers=(Layer(kind="uniform", background_material="dielectric",
                                        background_n=1.5, thickness_um=.1),),
                           incident_n=1, exit_n=1.5, order_budget=9, grid_size=16)
        result = linked_observable_sweep(model, "period_x_um",
            [{"path": "period_y_um", "scale": 1, "offset": 0}], .5, .7, 3,
            {"quantity": "T"})
        self.assertEqual(len(result["rows"]), 3)
        for row in result["rows"]:
            self.assertAlmostEqual(row["parameters"]["period_x_um"],
                                   row["parameters"]["period_y_um"])
            self.assertEqual(len(row["submitted_model_sha256"]), 64)

    def test_order_map_classifies_retained_orders(self):
        model = StackModel(layers=(Layer(kind="uniform", background_n=1.5),),
                           period_x_um=.4, period_y_um=.4, order_budget=9, grid_size=16)
        result = diffraction_order_map(model)
        zero = next(row for row in result["orders"] if row["m"] == row["n"] == 0)
        self.assertEqual(zero["ports"]["reflected"]["status"], "propagating")
        self.assertEqual(len(settings_fingerprint({"a": 1})), 64)


if __name__ == "__main__":
    unittest.main()
