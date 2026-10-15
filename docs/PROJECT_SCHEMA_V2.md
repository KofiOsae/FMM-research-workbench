# Project file schema version 2

Saved projects use `format: "fmm-research-workbench-project"` and `version: 2`.
The envelope keeps five concerns separate:

| Section | Meaning | Included in solver fingerprint? |
| --- | --- | --- |
| `solver_input` | Physical model, numerical settings, illumination direction, and active operation inputs | Yes |
| `workspace_state` | Selected workspace, tab, controls, and interface state | No |
| `result` | Completed result views and specialist-workspace records | No |
| `validation` | Numerical convergence and solver-specific evidence saved with the result | No |
| `environment` | Application release, Git commit, Python, platform, packages, and linear-algebra provider | No |

`solver_input_sha256` is the SHA-256 digest of recursively key-sorted,
compact UTF-8 JSON for `solver_input`. Opening a version-2 file recomputes this
digest and rejects a mismatch. The reference implementation is
`project_schema.py`.

Version-1 files remain readable. The browser migrates their layers, workspace,
controls, cached result views, and specialist data in memory. Since version 1
did not store a complete solver input or environment, the migrated file labels
those provenance gaps rather than inventing values. Saving it again produces a
version-2 file with a current fingerprint and environment record.

Changing column widths, focused/all-tools view, open cards, or other layout
preferences does not change the solver fingerprint. Changing a material,
geometry, illumination, numerical setting, or active calculation input does.
