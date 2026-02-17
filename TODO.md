# TODO - Think Like A Physicist (Quarto Book/Website)

## Status: Phase 1 in progress

---

## Decisions Made
- Repo root IS the Quarto project (restructured)
- Old attempts archived in `_archive/`
- IPhO note photos kept as reference only, not embedded in book
- Where diagrams are needed but not yet created, use a placeholder TODO callout
- **No labs or test-taking content** — per Casey's foreword, this covers thinking/physics only

---

### Questions for Casey
1. Error propagation — considered under lab stuff that we're not covering?

## Phase 1: Transcribe Casey's PDF Content → Quarto Chapters
- [x] **Foreword** (`index.qmd`) — Casey's foreword
- [x] **Chapter 1: How to Think Like a Physicist** (`thinking.qmd`)
  - [x] SI Units & dimensional analysis
  - [x] The 7 Ds problem-solving method
  - [x] Example Problem 1: Spring compression
  - [x] Example Problem 2: Piano on incline
  - [ ] Add placeholder diagrams for both example problems
- [x] **Chapter 2: Math Toolkit** (`math.qmd`)
  - [x] Differentiation (chain rule examples)
  - [x] Partial derivatives
  - [x] Integration (substitution, by parts)
  - [x] Kinematic equations
  - [x] Differential equations (1st order, SHM, separable)
  - [x] Coordinate systems (Cartesian, polar, cylindrical, spherical)
  - [x] Multivariate integration (areas, volumes, flux integrals)
- [x] ~~**Chapter 3: Labs & Error Analysis**~~ — Removed (out of scope per foreword)

## Phase 2: Transcribe IPhO Notes (169 photos → chapters)
- [x] Catalog all 169 photos by topic (see `images/IPhO_notes_catalog.md`)
- [ ] **Chapter 4: Mechanics** (`mechanics.qmd`) — forces, energy, momentum, rotation
- [ ] **Chapter 5: Thermodynamics** (`thermo.qmd`) — gas laws, entropy, heat engines
- [ ] **Chapter 6: Waves & Oscillations** (`waves.qmd`)
- [ ] **Chapter 7: Optics** (`optics.qmd`) — mirrors, lenses, diffraction, phasors
- [ ] **Chapter 8: Electrostatics** (`electrostatics.qmd`) — Coulomb, fields, Gauss's law
- [ ] **Chapter 9: Circuits & Magnetism** (`circuits.qmd`) — dipoles, induction, motors
- [ ] **Chapter 10: Modern Physics** (`relativity.qmd`) — relativity, light cones, Doppler

## Phase 3: Polish & Publish
- [ ] Create proper diagrams (TikZ or SVG) to replace placeholders
- [ ] Add practice problems at end of each chapter
- [ ] Write co-author preface (Jonathan + Casey)
- [ ] Review all math for correctness
- [ ] Test PDF output (ensure equations render properly)
- [ ] Deploy to GitHub Pages

---

## Book Outline

```
index.qmd              — Foreword (Casey)
part: "Foundations"
  thinking.qmd         — How to Think Like a Physicist (7 Ds, units, dimensions)
  math.qmd             — Math Toolkit (calculus, DEs, coordinates)
part: "Classical Physics"
  mechanics.qmd        — Mechanics
  thermo.qmd           — Thermodynamics
  waves.qmd            — Waves & Oscillations
part: "Electromagnetism & Optics"
  electrostatics.qmd   — Electrostatics
  circuits.qmd         — Circuits & Magnetism
  optics.qmd           — Optics
part: "Modern Physics"
  relativity.qmd       — Relativity
```

## Misc

1. Is there a way to have the calculations verifiable through something like sympy? It's nice to have equation manipulation through markdown/latex, but I am worried that it's easy to get a transcribed step that's wrong and hard to catch? Maybe not even in the published content but as a test process? Or do you have any similar suggestions?

---

## Open Questions
1. Any topics Casey wants to add beyond what's in the IPhO notes?
2. What order to tackle Phase 2 chapters?
3. Target audience refinement — "motivated professionals who know how to Google" per Casey's foreword, or broader?
