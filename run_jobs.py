"""Bounded local job queue and honest, cooperative calculation progress."""
from concurrent.futures import ThreadPoolExecutor
from contextvars import ContextVar
from threading import Lock
from time import monotonic, time
from uuid import uuid4
from pathlib import Path
import json

_current = ContextVar('calculation_job', default=None)
_lock = Lock()
_jobs = {}
_pool = ThreadPoolExecutor(max_workers=1, thread_name_prefix='solver')
_history_path = Path(__file__).with_name('.workbench_history.json')
_history_limit = 100


def _read_history():
    try:
        value = json.loads(_history_path.read_text(encoding='utf-8'))
        return value if isinstance(value, list) else []
    except (OSError, ValueError):
        return []


def _write_history():
    completed = []
    for key, job in _jobs.items():
        if job['status'] not in ('complete', 'failed', 'cancelled'):
            continue
        completed.append({
            'job_id': key, 'operation': job.get('operation', 'calculation'),
            'status': job['status'], 'created_at': job.get('created_at'),
            'finished_at': job.get('finished_at'),
            'elapsed_seconds': max(0., job.get('finished', monotonic())-
                                   job.get('started', job['created'])),
            'completed': job.get('completed', 0), 'total': job.get('total', 0),
            'stage': job.get('stage', ''),
            'error': (job.get('result') or {}).get('error')
        })
    prior = _read_history()
    known = {item.get('job_id') for item in completed}
    merged = completed + [item for item in prior if item.get('job_id') not in known]
    merged.sort(key=lambda item: item.get('finished_at') or item.get('created_at') or 0,
                reverse=True)
    try:
        _history_path.write_text(json.dumps(merged[:_history_limit], indent=2), encoding='utf-8')
    except OSError:
        pass


def progress(completed, total, stage='Calculating'):
    job = _current.get()
    if job is not None:
        if job.get('cancel_requested'):
            raise ValueError('Calculation cancelled by user')
        with _lock:
            job.update(completed=int(completed), total=int(total), stage=stage)


def submit(calculate, operation='calculation'):
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
        job = dict(status='queued', created=monotonic(), created_at=time(),
                   operation=str(operation), completed=0, total=0,
                   stage='Waiting for solver')
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
            job['finished_at'] = time()
            with _lock:
                _write_history()
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


def forget(key):
    """Release an acknowledged completed result while retaining history."""
    with _lock:
        job = _jobs.get(key)
        if job is None:
            return False
        if job['status'] in ('queued', 'running'):
            return False
        del _jobs[key]
        return True


def status_summary():
    """Small health payload that never exposes result arrays."""
    with _lock:
        states = [job['status'] for job in _jobs.values()]
    return {
        'queued': states.count('queued'),
        'running': states.count('running'),
        'retained_completed_results': sum(
            state in ('complete', 'failed', 'cancelled') for state in states),
        'maximum_active_jobs': 4,
        'solver_workers': 1,
    }


def history(limit=30):
    """Return recent durable job summaries, newest first."""
    limit = max(1, min(int(limit), _history_limit))
    active = []
    with _lock:
        for key, job in _jobs.items():
            if job['status'] in ('queued', 'running'):
                active.append({'job_id': key, 'operation': job.get('operation'),
                               'status': job['status'], 'created_at': job.get('created_at'),
                               'completed': job.get('completed', 0),
                               'total': job.get('total', 0), 'stage': job.get('stage', '')})
    active.sort(key=lambda item: item.get('created_at') or 0, reverse=True)
    completed = _read_history()
    completed.sort(key=lambda item: item.get('finished_at') or item.get('created_at') or 0,
                   reverse=True)
    return (active + completed)[:limit]
