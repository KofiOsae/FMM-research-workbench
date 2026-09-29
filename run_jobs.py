"""Bounded local job queue and honest, cooperative calculation progress."""
from concurrent.futures import ThreadPoolExecutor
from contextvars import ContextVar
from threading import Lock
from time import monotonic
from uuid import uuid4

_current = ContextVar('calculation_job', default=None)
_lock = Lock()
_jobs = {}
_pool = ThreadPoolExecutor(max_workers=1, thread_name_prefix='solver')


def progress(completed, total, stage='Calculating'):
    job = _current.get()
    if job is not None:
        if job.get('cancel_requested'):
            raise ValueError('Calculation cancelled by user')
        with _lock:
            job.update(completed=int(completed), total=int(total), stage=stage)


def submit(calculate):
    with _lock:
        # Keep result memory bounded; never evict an active calculation.
        for key in list(_jobs):
            if _jobs[key]['status'] in ('complete', 'failed', 'cancelled') and monotonic()-_jobs[key]['created'] > 1800:
                del _jobs[key]
        if sum(j['status'] in ('queued', 'running') for j in _jobs.values()) >= 4:
            raise ValueError('Four calculations are already queued; wait or cancel a queued run')
        if len(_jobs) >= 12:
            old = next((k for k,v in _jobs.items() if v['status'] in ('complete','failed','cancelled')), None)
            if old:
                del _jobs[old]
        key = uuid4().hex
        job = dict(status='queued', created=monotonic(), completed=0, total=0, stage='Waiting for solver')
        _jobs[key] = job

    def work():
        token = _current.set(job)
        try:
            if job.get('cancel_requested'):
                raise ValueError('Calculation cancelled by user')
            job.update(status='running', started=monotonic(), stage='Solving')
            result, status = calculate()
            job.update(result=result, response_status=status, status='complete' if status < 400 else 'failed')
        except Exception as exc:
            job.update(status='cancelled' if job.get('cancel_requested') else 'failed',
                       result={'error': str(exc)}, response_status=400)
        finally:
            if job.get('cancel_requested'):
                job.update(status='cancelled', result={'error':'Calculation cancelled by user'}, response_status=400)
            job['finished'] = monotonic()
            _current.reset(token)
    _pool.submit(work)
    return key


def snapshot(key, cancel=False):
    with _lock:
        job = _jobs.get(key)
        if job is None:
            return None
        if cancel and job['status'] in ('queued','running'):
            job['cancel_requested'] = True
        result = dict(job)
    elapsed = result.get('finished', monotonic())-result.get('started', result['created'])
    result['elapsed_seconds'] = max(0, elapsed)
    done, total = result['completed'], result['total']
    result['remaining_seconds'] = elapsed*(total-done)/done if done and total > done else None
    return {k:v for k,v in result.items() if k not in ('created','started','finished')}
