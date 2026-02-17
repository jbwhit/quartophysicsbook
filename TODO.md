# TODO - Think Like A Physicist (Quarto Book/Website)

## Status: Phase 2 in progress

---

## Decisions Made
- Repo root IS the Quarto project (restructured)
- Old attempts archived in `_archive/`
- IPhO note photos kept as reference only, not embedded in book
- Where diagrams are needed but not yet created, use a placeholder TODO callout
- **No labs or test-taking content** — per Casey's foreword, this covers thinking/physics only

---

## Phase 2: Transcribe IPhO Notes (169 photos → chapters)
- [x] Catalog all 169 photos by topic (see `images/IPhO_notes_catalog.md`)
- [x] **Chapter 4: Mechanics** (`mechanics.qmd`) — forces, energy, momentum, rotation
- [x] **Chapter 5: Thermodynamics** (`thermo.qmd`) — gas laws, entropy, heat engines
- [x] **Chapter 6: Waves & Oscillations** (`waves.qmd`)
- [x] **Chapter 7: Optics** (`optics.qmd`) — mirrors, lenses, diffraction, phasors
- [x] **Chapter 8: Electrostatics** (`electrostatics.qmd`) — Coulomb, fields, Gauss's law
- [x] **Chapter 9: Circuits & Magnetism** (`circuits.qmd`) — dipoles, induction, motors
- [ ] **Chapter 10: Modern Physics** (`relativity.qmd`) — relativity, light cones, Doppler

## Phase 3: Polish & Publish
- [ ] Create proper diagrams (TikZ or SVG) to replace 40 placeholders
- [ ] Add practice problems at end of each chapter
- [ ] Write co-author preface (Jonathan + Casey)
- [ ] Review all math for correctness
- [x] Test PDF output (renders cleanly; CI validates on every push)
- [x] Deploy to GitHub Pages (CI publishes HTML on push to `main`)

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

---

## Open Questions
1. Any topics Casey wants to add beyond what's in the IPhO notes?
2. Error propagation — considered under lab stuff that we're not covering?
