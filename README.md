# Think Like A Physicist

A terse undergraduate physics textbook by Casey Handmer and Jonathan Whitmore, built with [Quarto](https://quarto.org/).

Covers first-principles thinking and core physics — no lab sections, no test-taking content. Target audience: motivated professionals who know how to Google.

## Live site

Deployed to GitHub Pages on every push to `main`:
**https://jbwhit.github.io/quartophysicsbook/**

## Book structure

| Part | Chapter | File |
|------|---------|------|
| Foundations | How to Think Like a Physicist | `thinking.qmd` |
| | Math Toolkit | `math.qmd` |
| Classical Physics | Mechanics | `mechanics.qmd` |
| | Thermodynamics | `thermo.qmd` |
| | Waves & Oscillations | `waves.qmd` |
| Electromagnetism & Optics | Electrostatics | `electrostatics.qmd` |
| | Circuits & Magnetism | `circuits.qmd` |
| | Optics | `optics.qmd` |
| Modern Physics | Relativity | `relativity.qmd` |

## Getting started

### Prerequisites

- [Quarto](https://quarto.org/docs/get-started/) (1.4+)
- [uv](https://docs.astral.sh/uv/) (Python package manager)
- LaTeX (for PDF output): `quarto install tinytex`

### Preview locally

```bash
quarto preview
```

This starts a local server with live reload at `http://localhost:4200`.

### Render

```bash
# HTML site (output in _book/)
quarto render

# PDF only
quarto render --to pdf
```

### Run tests

Every boxed equation in the book has a corresponding [SymPy](https://www.sympy.org/) test that verifies the derivation symbolically.

```bash
uv run pytest tests/ -v
```

## Project conventions

- **Source material**: Casey's handwritten IPhO training note photos in `images/IPhO notes/` (reference only, not embedded in book). Catalog in `images/IPhO_notes_catalog.md`.
- **Worked problems** follow the **7 Ds** template: Diagram, Directions, Definitions, Diagnosis, Derivation, Determination, Dimensions (+ Substitution).
- **Diagram placeholders**: Where a diagram is needed but not yet created, a `{.callout-warning}` block describes what it should contain.
- **Tests**: One test file per chapter in `tests/test_<chapter>.py`. Run before committing.
- **Tracking**: `TODO.md` for current work, `DONE.md` for completed phases.

## CI/CD

Push to `main` triggers `.github/workflows/publish.yml`:
1. Renders PDF (catches any LaTeX errors)
2. Renders and publishes HTML to GitHub Pages

## License

All rights reserved. Contact the authors for permissions.
