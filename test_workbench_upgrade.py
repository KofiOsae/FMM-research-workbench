"""Regression checks for the training-driven workbench changes."""
import json
import tempfile
import time
import unittest
from dataclasses import replace
from unittest.mock import patch

import numpy as np
import run_jobs
from app import spectrum
from stack import StackModel, Layer
from tmm import solve_tmm
from scattering_maps import angle_wavelength_map
from bands import BandModel, full_zone_gaps


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

    def test_map_above_old_limit(self):
        with patch('scattering_maps.solve_stack',return_value={'R':.04}), patch('scattering_maps.stack_convergence',return_value={'converged':True}):
            result=angle_wavelength_map(self.model(),.8,1.2,43,-5,5,3)
        self.assertEqual(np.shape(result['values']),(3,43))

    def test_full_zone_gap_screen(self):
        result=full_zone_gaps(BandModel(background_n=2,inclusion_n=2,
            fourier_order=1,grid_size=32,points_per_segment=3,bands=3),5)
        self.assertEqual(result['grid_points_per_axis'],5)
        self.assertEqual(result['complete_sampled_gaps'],[])

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
            with urlopen(base+'/workbench.js') as response:
                script=response.read()
                self.assertIn(b'Run from here',script)
                self.assertIn(b'Straight line',script)
                self.assertIn(b'materialPreviewColors',script)
                self.assertIn(b'periodic pixel-run screen',script)
                self.assertIn(b'Export SVG',script)
                self.assertIn(b'maskBoundarySegments',script)
                self.assertIn(b'extruded through the entered layer thickness',script)
                self.assertIn(b'Geometry in-plane rotation',script)
                self.assertIn(b'rotateMask(angleDeg)',script)
                self.assertIn(b'Add transformed copy',script)
                self.assertIn(b'transformedPixels',script)
            with urlopen(base+'/METASURFACE_GUIDE.html') as response:
                self.assertIn(b'Flat-lens design', response.read())
            with urlopen(base+'/') as response:
                page=response.read()
                self.assertIn(b'Tilt \xe2\x86\x91',page)
                self.assertIn('k∥ / k₀'.encode(),page)
                self.assertIn(b'Previous successful result',page)
                self.assertIn(b'run-diagnostic',page)
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
