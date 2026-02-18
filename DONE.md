# DONE - Think Like A Physicist

Completed work, moved here from TODO.md to keep it clean.

---

## Phase 0: Project Setup & Structure (COMPLETE)
- [x] Restructure repo root as Quarto project
- [x] Archive old attempts (`firstdraft/`, `old-qph/`, `physicsthinking/`, `repo/`, `quartophysicsbook/`) into `_archive/`
- [x] Set up `_quarto.yml` with dual output (website + PDF book)
- [x] Set up GitHub Actions deploy workflow
- [x] Create chapter structure with parts
- [x] Initialize git repo at root, initial commit, push to GitHub
- [x] Verify `quarto preview` works locally

## Phase 1: Transcribe Casey's PDF Content (COMPLETE)
- [x] **Foreword** (`index.qmd`) — Casey's foreword
- [x] **Chapter 1: How to Think Like a Physicist** (`thinking.qmd`)
  - [x] SI Units & dimensional analysis
  - [x] The 7 Ds problem-solving method
  - [x] Example Problem 1: Spring compression
  - [x] Example Problem 2: Piano on incline
  - [x] Placeholder diagrams for both example problems
- [x] **Chapter 2: Math Toolkit** (`math.qmd`)
  - [x] Differentiation (chain rule examples)
  - [x] Partial derivatives
  - [x] Integration (substitution, by parts)
  - [x] Kinematic equations
  - [x] Differential equations (1st order, SHM, separable)
  - [x] Coordinate systems (Cartesian, polar, cylindrical, spherical)
  - [x] Multivariate integration (areas, volumes, flux integrals)
- [x] ~~**Chapter 3: Labs & Error Analysis**~~ — Removed (out of scope per foreword)

## Phase 2: Transcribe IPhO Notes — 169 photos (COMPLETE)
- [x] Catalog all 169 photos by topic (see `images/IPhO_notes_catalog.md`)
- [x] **Chapter 4: Mechanics** (`mechanics.qmd`)
- [x] **Chapter 5: Thermodynamics** (`thermo.qmd`)
- [x] **Chapter 6: Waves & Oscillations** (`waves.qmd`)
- [x] **Chapter 7: Optics** (`optics.qmd`)
- [x] **Chapter 8: Electrostatics** (`electrostatics.qmd`)
- [x] **Chapter 9: Circuits & Magnetism** (`circuits.qmd`)
- [x] **Chapter 10: Modern Physics & Relativity** (`relativity.qmd`)

## Infrastructure (COMPLETE)
- [x] Test PDF output (renders cleanly; CI validates on every push)
- [x] Deploy to GitHub Pages (CI publishes HTML on push to `main`)
- [x] Install `pandoc-ext/diagram` extension for TikZ support
- [x] Set up CI with TinyTeX + Inkscape for diagram rendering
- [x] Add "View source" and "Report an issue" links to website

## Decisions Made
- Repo root IS the Quarto project (restructured)
- Old attempts archived in `_archive/`
- IPhO note photos kept as reference only, not embedded in book
- Where diagrams are needed but not yet created, use a placeholder TODO callout
- **No labs or test-taking content** — per Casey's foreword, this covers thinking/physics only
- AI-generated TikZ diagrams get a visible "Needs Review" callout until human sign-off
