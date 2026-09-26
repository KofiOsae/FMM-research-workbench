# Contributing

Thank you for improving FMM Research Workbench. Changes that affect scientific results must include evidence appropriate to the claim.

## Before coding

1. Search existing issues.
2. Open an issue describing the problem, expected result, literature or analytic reference, and a minimal exported project file.
3. Keep one pull request focused on one change.

## Local setup

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m unittest test_validation.py
.\.venv\Scripts\python.exe app.py
```

Open `http://127.0.0.1:8765/`.

## Scientific changes

Include at least one of these forms of evidence:

- an analytic limiting case;
- comparison with an independent solver;
- a digitized, correctly attributed literature benchmark;
- convergence under Fourier order, geometry grid, field grid, or domain refinement.

State the observable and convention precisely: total or specular power, polarization basis, angle convention, normalization, material specimen, and wavelength units. Visual similarity alone is not numerical validation.

## Pull requests

- Add or update a meaningful test.
- Run the complete test suite.
- Update the tutorial and validation notes when capabilities or limits change.
- Record new material provenance, wavelength range, temperature, specimen, and whether absorption is included.
- Do not commit measurements, proprietary data, credentials, logs, virtual environments, or local `user_materials.json` files.

By submitting a contribution, you agree that it is distributed under the repository license.
