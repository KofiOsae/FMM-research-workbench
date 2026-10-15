"""Versioned project envelope, migration, and solver-input identity.

Only ``solver_input`` contributes to the scientific fingerprint. Workspace
layout, cached results, validation notes, and runtime provenance remain
separate so rearranging the interface cannot change model identity.
"""

from __future__ import annotations

import hashlib
import json
from copy import deepcopy
from typing import Any

PROJECT_FORMAT = "fmm-research-workbench-project"
CURRENT_VERSION = 2


def canonical_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True,
                      separators=(",", ":"), allow_nan=False)


def solver_input_fingerprint(solver_input: dict) -> str:
    if not isinstance(solver_input, dict):
        raise ValueError("solver_input must be an object")
    return hashlib.sha256(canonical_json(solver_input).encode("utf-8")).hexdigest()


def migrate_project(project: dict) -> dict:
    if not isinstance(project, dict) or project.get("format") != PROJECT_FORMAT:
        raise ValueError("This is not an FMM Research Workbench project")
    version = project.get("version", 1)
    if version == CURRENT_VERSION:
        return deepcopy(project)
    if version != 1 or not isinstance(project.get("layers"), list):
        raise ValueError(f"Unsupported project version: {version}")
    return {
        "format": PROJECT_FORMAT, "version": CURRENT_VERSION,
        "saved_at": project.get("savedAt"),
        "solver_input": {
            "model": None, "physical_layers": deepcopy(project["layers"]),
            "active_operation": {"kind": project.get("tab", "spectrum"), "inputs": {}},
        },
        "solver_input_sha256": None,
        "workspace_state": {
            "workspace": project.get("workspace", "scattering"),
            "tab": project.get("tab", "spectrum"),
            "controls": deepcopy(project.get("inputs", {})),
        },
        "result": {"views": deepcopy(project.get("resultViews", {})),
                   "specialists": deepcopy(project.get("extras", {}))},
        "validation": {"migration": "Legacy version 1 project; no solver fingerprint was stored."},
        "environment": {"migration": "Runtime provenance was unavailable in version 1."},
        "migration": {"from_version": 1},
    }


def validate_project(project: dict) -> dict:
    normalized = migrate_project(project)
    required = ("solver_input", "workspace_state", "result", "validation", "environment")
    missing = [name for name in required if not isinstance(normalized.get(name), dict)]
    if missing:
        raise ValueError("Project version 2 is missing object sections: " + ", ".join(missing))
    if not isinstance(normalized["solver_input"].get("physical_layers"), list):
        raise ValueError("solver_input.physical_layers must be an array")
    actual = solver_input_fingerprint(normalized["solver_input"])
    stored = normalized.get("solver_input_sha256")
    if stored is not None and stored != actual:
        raise ValueError("Project solver-input fingerprint does not match its contents")
    normalized["solver_input_sha256"] = actual
    return normalized
