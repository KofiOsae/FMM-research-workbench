"""Regression checks for the training-driven workbench changes."""
import json
import tempfile
import time
import unittest
from pathlib import Path
from dataclasses import replace
from unittest.mock import patch

import numpy as np
import run_jobs
import usage_stats
from app import spectrum, _spectrum_diagnostics
from stack import StackModel, Layer
from tmm import solve_tmm
from scattering_maps import angle_wavelength_map
from bands import BandModel, full_zone_gaps
from sweep import diffraction_order_sweep
from finite_grating import FiniteGratingModel, solve_finite_grating
from finite_grating_bidirectional import reciprocity_certificate


class UpgradeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._history_directory = tempfile.TemporaryDirectory()
        cls._original_history_path = run_jobs._history_path
        run_jobs._history_path = run_jobs.Path(cls._history_directory.name)/'history.json'

    @classmethod
    def tearDownClass(cls):
        run_jobs._history_path = cls._original_history_path
        cls._history_directory.cleanup()

    def model(self):
        return StackModel(incident_n=1, exit_n=1.5, wavelength_um=1,
                          layers=(Layer(kind='uniform', thickness_um=.125,
                                        background_material='dielectric', background_n=2),))

    def test_thin_film_handbook_reference(self):
        result = solve_tmm(self.model())
        self.assertAlmostEqual(result['R'], .206611570247934, places=12)
        self.assertAlmostEqual(result['R']+result['T'], 1, places=12)
        self.assertAlmostEqual(abs(result['r_phase_deg']), 180, places=10)
        self.assertAlmostEqual(result['t_phase_deg'], 90, places=10)

    def test_unpolarized_has_no_average_phase(self):
        result=solve_tmm(replace(self.model(), polarization='unpolarized'))
        self.assertIsNone(result['r_phase_deg'])
        self.assertIn('s',result['polarization_channels'])

    def test_high_resolution_spectrum(self):
        self.assertEqual(len(spectrum(self.model(), .8, 1.2, 151)['rows']),151)

    def test_spectrum_failures_report_cause_and_action(self):
        patterned=replace(self.model(),layers=(Layer(kind='stripe',thickness_um=.125,
                            background_material='air',feature_material='dielectric',feature_n=2),))
        with patch('app.stack_convergence', side_effect=ValueError('An order is at grazing cutoff; adjust wavelength or angle slightly')):
            result=spectrum(patterned, .8, .82, 3)
        diagnostics=result['diagnostics']
        self.assertEqual(diagnostics['failed_points'],3)
        self.assertEqual(diagnostics['computed_points'],0)
        self.assertIn('grazing cutoff',diagnostics['causes'][0]['message'])
        self.assertIn('Move the wavelength range',diagnostics['causes'][0]['suggestion'])
        material=_spectrum_diagnostics([{'wavelength_um':5.0,'status':'Green silicon data cover 0.25–1.45 µm'}])
        self.assertIn('material-data range',material['causes'][0]['suggestion'])

    def test_map_above_old_limit(self):
        with patch('scattering_maps.solve_stack',return_value={'R':.04}), patch('scattering_maps.stack_convergence',return_value={'converged':True}):
            result=angle_wavelength_map(self.model(),.8,1.2,43,-5,5,3)
        self.assertEqual(np.shape(result['values']),(3,43))

    def test_full_zone_gap_screen(self):
        result=full_zone_gaps(BandModel(background_n=2,inclusion_n=2,
            fourier_order=1,grid_size=32,points_per_segment=3,bands=3),5)
        self.assertEqual(result['grid_points_per_axis'],5)
        self.assertEqual(result['complete_sampled_gaps'],[])

    def test_diffraction_order_sweep_opens_first_orders(self):
        grating=StackModel(wavelength_um=.633,incident_n=1,exit_n=1,
            period_x_um=.55,period_y_um=.3,order_budget=13,grid_size=24,
            layers=(Layer(kind='stripe',thickness_um=.2,
                background_material='air',feature_material='dielectric',
                feature_n=2,fill_x=.5,fill_y=.5),))
        result=diffraction_order_sweep(grating,'period_x_um',.55,1.05,3)
        self.assertEqual(result['solved_points'],3)
        self.assertIn({'m':-1,'n':0},result['orders'])
        self.assertIn({'m':1,'n':0},result['orders'])
        self.assertEqual({(o['m'],o['n']) for o in result['rows'][0]['orders']},{(0,0)})
        last=result['rows'][-1]
        self.assertAlmostEqual(sum(o['R'] for o in last['orders']),last['R'],places=10)
        self.assertAlmostEqual(sum(o['T'] for o in last['orders']),last['T'],places=10)

    def test_finite_grating_has_two_independent_reciprocal_sources(self):
        model=FiniteGratingModel(mesh_um=.1,periods=4,left_padding_um=1,
            right_padding_um=1.2,top_padding_um=1.2,substrate_depth_um=1.2,
            absorber_um=.3,target_waist_um=1)
        incoming=solve_finite_grating(replace(model,excitation='fiber'))
        self.assertEqual(incoming['formulation']['excitation'],
                         'gaussian_port_to_grating_to_waveguide')
        self.assertLess(incoming['formulation']['reference_reproduction_relative_error'],1e-10)
        self.assertIn('incident_Ey2',incoming['field'])
        self.assertIn('scattered_Ey2',incoming['field'])
        self.assertIn('incident',incoming['angular_spectrum'])
        self.assertIn('annotations',incoming['field'])
        certificate=reciprocity_certificate(model)
        self.assertLess(certificate['relative_magnitude_error'],.1)
        self.assertGreater(certificate['forward_complex_coefficient']['magnitude'],0)
        self.assertGreater(certificate['reverse_complex_coefficient']['magnitude'],0)

    def wait(self,key):
        deadline=time.monotonic()+3
        while time.monotonic()<deadline:
            value=run_jobs.snapshot(key)
            if value['status'] in ('complete','failed','cancelled'):
                return value
            time.sleep(.01)
        self.fail('Job failed to finish')

    def test_job_completion_and_progress(self):
        def work():
            run_jobs.progress(2,3,'Test')
            return {'value':42},200
        result=self.wait(run_jobs.submit(work))
        self.assertEqual(result['result']['value'],42)
        self.assertEqual(result['completed'],2)
        json.dumps(result,allow_nan=False)
        self.assertTrue(any(item['status']=='complete' for item in run_jobs.history()))

    def test_completed_job_result_can_be_released(self):
        key=run_jobs.submit(lambda: ({'large_result':[1,2,3]},200))
        self.assertEqual(self.wait(key)['status'],'complete')
        self.assertTrue(run_jobs.forget(key))
        self.assertIsNone(run_jobs.snapshot(key))
        self.assertTrue(any(item['job_id']==key for item in run_jobs.history()))

    def test_anonymous_usage_statistics_are_aggregate_and_deduplicated(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(
                usage_stats, '_path', Path(directory) / 'usage.json'):
            first = usage_stats.record_visit('random-browser-12345', True)
            second = usage_stats.record_visit('random-browser-12345', False)
            usage_stats.record_calculation('spectrum', 'started')
            usage_stats.record_calculation('spectrum', 'completed')
            usage_stats.record_calculation('bands', 'started')
            usage_stats.record_calculation('bands', 'failed')
            result = usage_stats.summary()
        self.assertTrue(first['new_anonymous_visitor'])
        self.assertFalse(second['new_anonymous_visitor'])
        self.assertEqual(result['totals']['anonymous_visitors'], 1)
        self.assertEqual(result['totals']['page_views'], 2)
        self.assertEqual(result['totals']['sessions'], 1)
        self.assertEqual(result['totals']['calculations_completed'], 1)
        self.assertEqual(result['totals']['calculations_failed'], 1)
        self.assertNotIn('random-browser-12345', json.dumps(result))

    def test_job_failure(self):
        def work():
            raise ValueError('Known failure')
        result=self.wait(run_jobs.submit(work))
        self.assertEqual(result['status'],'failed')
        self.assertEqual(result['result']['error'],'Known failure')

    def test_job_cancellation(self):
        from threading import Event
        release=Event()
        def work():
            release.wait(2)
            run_jobs.progress(1,2)
            return {},200
        key=run_jobs.submit(work)
        run_jobs.snapshot(key,cancel=True)
        release.set()
        self.assertEqual(self.wait(key)['status'],'cancelled')

    def test_http_queue_and_assets(self):
        from http.server import ThreadingHTTPServer
        from threading import Thread
        from urllib.request import Request, urlopen
        from dataclasses import asdict
        from app import Handler
        class QuietHandler(Handler):
            def log_message(self,*args):
                pass
        server=ThreadingHTTPServer(('127.0.0.1',0),QuietHandler)
        thread=Thread(target=server.serve_forever,daemon=True);thread.start()
        base=f'http://127.0.0.1:{server.server_port}'
        try:
            with urlopen(base+'/health') as response:
                health=json.loads(response.read())
                self.assertEqual(health['status'],'ok')
                self.assertEqual(health['version'],'bidirectional-grating-2026-10-10')
                self.assertEqual(health['queue']['solver_workers'],1)
                self.assertEqual(health['limits']['finite_grid_cells'],300000)
            with urlopen(base+'/api/usage-stats') as response:
                statistics=json.loads(response.read())
                self.assertIn('anonymous_visitors',statistics['totals'])
                self.assertIn('No IP address',statistics['privacy'])
            with urlopen(base+'/workbench.js') as response:
                script=response.read()
                self.assertIn(b'Run from here',script)
                self.assertIn(b'Straight line',script)
                self.assertIn(b'materialPreviewColors',script)
                self.assertIn(b"addEventListener('input',updateLegend)",script)
                self.assertIn(b'drawImportedMasks();updateLegend()',script)
                self.assertIn(b'periodic pixel-run screen',script)
                self.assertIn(b'Export SVG',script)
                self.assertIn(b'maskBoundarySegments',script)
                self.assertIn(b'extruded through the entered layer thickness',script)
                self.assertIn(b'Geometry in-plane rotation',script)
                self.assertIn(b'rotateMask(angleDeg)',script)
                self.assertIn(b'Add transformed copy',script)
                self.assertIn(b'transformedPixels',script)
                self.assertIn(b'Diffraction-order efficiency sweep',script)
                self.assertIn(b'runDiffractionSweep',script)
                self.assertIn(b"operation==='spectrum'",script)
                self.assertIn(b'collectWorkbenchProjectExtras',script)
                self.assertIn(b'Gaussian or imported fiber illumination',script)
                self.assertIn(b'parsePortCsv',script)
                self.assertIn(b'forward_efficiency',script)
                self.assertIn(b'Uniform-port source',script)
                self.assertIn(b'Validate mesh',script)
                self.assertIn(b'overlap_density_normalized',script)
                self.assertIn('Publication · 2×10⁻⁴'.encode(),script)
                self.assertIn(b'validationProfile',script)
                self.assertIn(b'portValidationProfile',script)
                self.assertIn(b'Upper/lower coupled-branch composition',script)
                self.assertIn(b'useTrackedBranches',script)
                self.assertIn(b'Reusable optical observables',script)
                self.assertIn(b'optimizerOpticalConstraints',script)
                self.assertIn(b'Design-result dashboard',script)
                self.assertIn(b'Finite grating',script)
                self.assertIn(b'finite_grating_coupler',script)
                self.assertIn(b'runFiniteGrating',script)
                self.assertIn(b'finite_grating_benchmark',script)
                self.assertIn(b'finite_grating_reciprocity',script)
                self.assertIn(b'Gaussian port',script)
                self.assertIn(b'incident_Ey2',script)
                self.assertIn(b'finiteAngularPlot',script)
                self.assertIn(b'Power fraction density',script)
                self.assertIn(b'runFgSpectrum',script)
                self.assertIn(b'runFgOptimize',script)
                self.assertIn(b'runFgTolerance',script)
                self.assertIn(b'finite_grating_device',script)
                self.assertIn(b'finiteGratingSchematic',script)
                self.assertIn(b'Parameter definitions, supported coupling, and limits',script)
                self.assertIn(b'How to read this result',script)
                self.assertIn(b'fgFiniteBox',script)
                self.assertIn(b'fgDiscretization',script)
                self.assertIn(b'fgFieldQuantity',script)
                self.assertIn(b'fgGridEstimate',script)
                self.assertIn(b'updateFiniteGridEstimate',script)
                self.assertIn(b'normalized_Ey_real',script)
                self.assertIn(b'Usage statistics',script)
                self.assertIn(b'Anonymous browsers',script)
                self.assertIn(b'/api/usage/visit',script)
                self.assertIn(b'Workspace layout',script)
                self.assertIn(b'workspaceLayout',script)
            with urlopen(base+'/workbench.css') as response:
                stylesheet=response.read()
                self.assertIn(b'orientation:landscape',stylesheet)
                self.assertIn(b'workspace-layout="split"',stylesheet)
                self.assertIn(b'data-workspace="usage"',stylesheet)
                self.assertIn(b'minmax(300px,340px)',stylesheet)
            with urlopen(base+'/METASURFACE_GUIDE.html') as response:
                guide=response.read()
                self.assertIn(b'Flat-lens design',guide)
                self.assertIn(b'Geometry in-plane rotation',guide)
                self.assertIn(b'Add a transformed copy in this layer',guide)
            with urlopen(base+'/FIRST_STEPS.html') as response:
                guide=response.read()
                self.assertIn(b'Diffraction-order efficiency sweep',guide)
                self.assertIn(b'Rayleigh threshold',guide)
                self.assertIn(b'Uniform-port projection of a Gaussian or imported fiber field',guide)
                self.assertIn(b'Uniform-port projection',guide)
                self.assertIn(b'Fit upper and lower coupled branches',guide)
                self.assertIn(b'Optimize a 50:50 two-order splitter',guide)
                self.assertIn(b'From a physical claim to a publication package',guide)
                self.assertIn(b'Publication uses 2&times;10<sup>-4</sup>',guide)
                self.assertIn(b'Calculate grating coupling in both reciprocal directions',guide)
                self.assertIn(b'Research map for each capability',guide)
                self.assertIn(b'Marchetti et al.',guide)
            with urlopen(base+'/') as response:
                page=response.read()
                self.assertIn(b'Tilt \xe2\x86\x91',page)
                self.assertIn('k∥ / k₀'.encode(),page)
                self.assertIn(b'Previous successful result',page)
                self.assertIn(b'run-diagnostic',page)
                self.assertIn(b'appendSpectrumDiagnostics',page)
                self.assertIn(b'Primary solver message',page)
                self.assertIn(b'restoreWorkbenchProjectExtras',page)
                self.assertIn('Reflected order Rₘₙ'.encode(),page)
                self.assertIn(b'toleranceOperatingLambda',page)
                self.assertIn(b'max="201" value="37"',page)
                self.assertIn(b'max="512" value="64"',page)
            with urlopen(base+'/workbench.js') as response:
                enhancements=response.read()
                self.assertIn(b'Validation report and diffraction light cones',enhancements)
                self.assertIn(b'Linked custom-observable sweep',enhancements)
                self.assertIn(b'Find a tool in this workspace',enhancements)
            body=json.dumps({'model':asdict(self.model()),'start':.8,'stop':1.2,'points':151}).encode()
            request=Request(base+'/api/tmm',data=body,headers={'Content-Type':'application/json','X-Workbench-Job':'1'})
            with urlopen(request) as response:
                self.assertEqual(response.status,202)
                key=json.load(response)['job_id']
            result=self.wait(key)
            self.assertEqual(result['status'],'complete')
            with urlopen(base+'/api/jobs/'+key) as response:
                self.assertEqual(len(json.load(response)['result']['rows']),151)
            with urlopen(base+'/api/job-history') as response:
                jobs=json.load(response)['jobs']
                self.assertTrue(any(item['job_id']==key and item['operation']=='tmm' for item in jobs))
        finally:
            server.shutdown();server.server_close();thread.join()

if __name__=='__main__':
    unittest.main()
