"""Privacy-preserving aggregate usage counters for the Workbench.

The browser creates a random identifier. Only a short one-way hash is stored;
IP addresses, user agents, inputs, results, and uploaded data are never read.
"""
from __future__ import annotations

from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from threading import Lock
import json
import os


_lock = Lock()
_path = Path(os.environ.get("USAGE_STATS_PATH",
                            Path(__file__).with_name(".usage_stats.json")))
_excluded_operations = {"materials", "material_import", "figure", "research_report"}


def _now():
    return datetime.now(timezone.utc)


def _blank():
    now = _now().isoformat()
    return {
        "schema_version": 1,
        "started_utc": now,
        "updated_utc": now,
        "anonymous_visitors": [],
        "totals": {"page_views": 0, "sessions": 0,
                   "calculations_started": 0, "calculations_completed": 0,
                   "calculations_failed": 0},
        "operations": {},
        "days": {},
    }


def _read():
    try:
        value = json.loads(_path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else _blank()
    except (OSError, ValueError):
        return _blank()


def _write(value):
    value["updated_utc"] = _now().isoformat()
    try:
        _path.parent.mkdir(parents=True, exist_ok=True)
        temporary = _path.with_suffix(_path.suffix + ".tmp")
        temporary.write_text(json.dumps(value, indent=2, sort_keys=True), encoding="utf-8")
        temporary.replace(_path)
        return True
    except OSError:
        # Usage measurement must never prevent a scientific calculation.
        return False


def _token(value: str) -> str:
    text = str(value or "")[:200]
    return sha256(("fmm-workbench-anonymous-v1:" + text).encode()).hexdigest()[:24]


def _day(value, date):
    days = value.setdefault("days", {})
    row = days.setdefault(date, {"page_views": 0, "sessions": 0,
                                 "anonymous_visitors": [],
                                 "calculations_completed": 0,
                                 "calculations_failed": 0})
    # Daily visitor hashes are retained for at most 90 days.
    for old in sorted(days)[:-90]:
        del days[old]
    return row


def record_visit(visitor_id: str, new_session: bool = True):
    if not isinstance(visitor_id, str) or not 8 <= len(visitor_id) <= 200:
        raise ValueError("Anonymous visitor identifier is invalid")
    token = _token(visitor_id)
    today = _now().date().isoformat()
    with _lock:
        value = _read()
        visitors = value.setdefault("anonymous_visitors", [])
        is_new_visitor = token not in visitors
        if is_new_visitor:
            visitors.append(token)
        totals = value.setdefault("totals", {})
        totals["page_views"] = int(totals.get("page_views", 0)) + 1
        if new_session:
            totals["sessions"] = int(totals.get("sessions", 0)) + 1
        row = _day(value, today)
        row["page_views"] = int(row.get("page_views", 0)) + 1
        if new_session:
            row["sessions"] = int(row.get("sessions", 0)) + 1
        daily = row.setdefault("anonymous_visitors", [])
        if token not in daily:
            daily.append(token)
        _write(value)
    return {"recorded": True, "new_anonymous_visitor": is_new_visitor}


def tracks_calculation(operation: str) -> bool:
    return bool(operation and "/" not in operation and
                operation not in _excluded_operations)


def record_calculation(operation: str, status: str):
    if not tracks_calculation(operation) or status not in {"started", "completed", "failed"}:
        return
    today = _now().date().isoformat()
    with _lock:
        value = _read()
        totals = value.setdefault("totals", {})
        key = "calculations_" + status
        totals[key] = int(totals.get(key, 0)) + 1
        operations = value.setdefault("operations", {})
        item = operations.setdefault(operation, {"started": 0, "completed": 0, "failed": 0})
        item[status] = int(item.get(status, 0)) + 1
        if status in {"completed", "failed"}:
            row = _day(value, today)
            daily_key = "calculations_" + status
            row[daily_key] = int(row.get(daily_key, 0)) + 1
        _write(value)


def summary(days: int = 30):
    days = max(1, min(int(days), 90))
    with _lock:
        value = _read()
    totals = dict(value.get("totals", {}))
    totals["anonymous_visitors"] = len(value.get("anonymous_visitors", []))
    daily = []
    for date, row in sorted(value.get("days", {}).items())[-days:]:
        daily.append({"date": date,
                      "anonymous_visitors": len(row.get("anonymous_visitors", [])),
                      "sessions": int(row.get("sessions", 0)),
                      "page_views": int(row.get("page_views", 0)),
                      "calculations_completed": int(row.get("calculations_completed", 0)),
                      "calculations_failed": int(row.get("calculations_failed", 0))})
    operations = [{"operation": name, **counts}
                  for name, counts in value.get("operations", {}).items()]
    operations.sort(key=lambda row: row.get("completed", 0), reverse=True)
    return {"schema_version": value.get("schema_version", 1),
            "started_utc": value.get("started_utc"),
            "updated_utc": value.get("updated_utc"),
            "totals": totals, "daily": daily,
            "operations": operations,
            "privacy": "Anonymous random browser identifiers are hashed. No IP address, user agent, uploaded data, simulation inputs, or results are stored.",
            "persistence": "Statistics persist only as long as the configured USAGE_STATS_PATH storage persists. Configure a durable mounted path for lifetime production totals."}
