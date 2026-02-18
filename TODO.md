# TODO - Think Like A Physicist

## Current State (for Casey sync, Feb 2026)

### Book Structure

| Part | Chapter | File | Lines | Boxed Results | Tests | TikZ Diagrams | TODO Diagrams |
|------|---------|------|-------|---------------|-------|---------------|---------------|
| Foundations | How to Think Like a Physicist | `thinking.qmd` | 177 | 2 | 6 | 0 | 3 |
| Foundations | Math Toolkit | `math.qmd` | 159 | 4 | 5 | 0 | 0 |
| Classical | Mechanics | `mechanics.qmd` | 877 | 20 | 21 | 0 | 8 |
| Classical | Thermodynamics | `thermo.qmd` | 486 | 11 | 14 | 1 | 6 |
| Classical | Waves & Oscillations | `waves.qmd` | 470 | 16 | 21 | 1 | 7 |
| E&M/Optics | Electrostatics | `electrostatics.qmd` | 279 | 15 | 25 | 0 | 2 |
| E&M/Optics | Circuits & Magnetism | `circuits.qmd` | 506 | 34 | 36 | 0 | 5 |
| E&M/Optics | Optics | `optics.qmd` | 511 | 16 | 29 | 0 | 11 |
| Modern | Modern Physics & Relativity | `relativity.qmd` | 443 | 18 | 35 | 1 | 3 |
| **Totals** | | | **~3,900** | **136** | **192** | **3** | **45** |

### Content Summary Per Chapter

- **Thinking**: SI units, dimensional analysis, 7 Ds method, 2 worked examples (spring compression, piano on incline)
- **Math**: Chain rule, partial derivatives, integration (substitution, by parts), kinematic equations, DEs (1st order, SHM), coordinate systems, multivariate integration
- **Mechanics**: Energy & conservation, Newton's laws, FBDs, friction, momentum, collisions (elastic/inelastic), rotational mechanics (torque, moment of inertia, rolling, angular momentum, precession), gravity & orbits (Kepler's laws, escape velocity), fluids (pressure, Archimedes, Bernoulli, Poiseuille). 2 full 7-Ds worked problems (Atwood machine, satellite orbit)
- **Thermo**: Temperature, heat transfer (conduction/radiation/convection), thermal expansion, phase changes, ideal gas law, kinetic theory, equipartition, thermodynamic processes (isochoric/isobaric/isothermal/adiabatic), heat engines, Carnot cycle, entropy, phase diagrams
- **Waves**: SHM, pendulum, traveling waves, energy/power in waves, superposition, standing waves, resonance, beats, Doppler effect. 2 full 7-Ds worked problems (guitar string, beat frequency)
- **Electrostatics**: Coulomb's law, electric field, Gauss's law (point charge, line charge, conducting sphere), electric potential, capacitors, dielectrics, electric dipoles
- **Circuits**: Ohm's law, Kirchhoff's laws, resistors/capacitors/inductors (series & parallel), RC/RL/LC transients, magnetostatics (Biot-Savart, Ampere's law, solenoids), electromagnetic induction (Faraday, Lenz), AC circuits (impedance, resonance, power, RMS), Maxwell's equations, EM waves, Poynting vector
- **Optics**: Reflection, Snell's law, total internal reflection, mirrors (concave/convex), thin lenses, lensmaker's equation, the eye, telescopes, double slit interference, diffraction gratings, thin film interference, Michelson interferometer, single slit diffraction, Airy disk, Rayleigh criterion, polarization (Malus's law, Brewster's angle, birefringence)
- **Relativity**: Bohr atom, photoelectric effect, wave-particle duality, de Broglie wavelength, uncertainty principle, special relativity (Lorentz transforms, length contraction, time dilation, velocity addition), spacetime diagrams, relativistic dynamics (four-vectors, E=mc^2, energy-momentum relation), Compton scattering, aberration, relativistic Doppler

### What's Done

- All 169 IPhO training note photos transcribed and catalogued
- All content from Casey's "Think Like A Physicist" PDF transcribed (Chapters 1-2)
- 136 boxed results with 192 passing sympy verification tests
- 3 AI-generated TikZ diagrams (pendulum, heat engine, light clock) — awaiting human review
- Physics correctness pass complete (5 errors fixed: PSI conversion, gravitational force mass, optics phase change rule, RC circuit equation, Doppler angle description)
- Consistency pass complete (vector notation unified to `\vec{}`, SI units throughout, `\delta`/`d`/`\partial` usage verified)
- CI/CD pipeline: GitHub Actions builds HTML + PDF and deploys to GitHub Pages on push to main
- Website live with dual light/dark theme and PDF download link in sidebar
- Pre-push git hook renders PDF locally before each push

### What's Left

- 45 TODO diagram placeholders need TikZ implementations
- 3 existing TikZ diagrams need human review
- Only 4 chapters have worked problems (thinking: 2, mechanics: 2, waves: 2); remaining 5 chapters have none
- No preface yet
- No cross-references between chapters

### Key Questions for Casey

1. **Topic coverage**: Are there topics Casey wants to add beyond what's in the IPhO notes? (e.g., error propagation was excluded as "lab stuff" — agree?)
2. **Depth**: Is the level of detail right? Some chapters are dense (circuits: 506 lines, 34 boxed results) while others are lighter (electrostatics: 279 lines)
3. **Tone**: The transcribed IPhO notes are terse and formula-heavy. Casey's PDF chapters (thinking, math) are more conversational. Should we unify the tone?
4. **Worked problems**: Currently only in thinking/mechanics/waves. Should every chapter get 7-Ds practice problems?
5. **Modern physics scope**: The relativity chapter covers Bohr atom, photoelectric effect, QM intro, and SR. Is this the right mix for the book's level?
6. **Diagrams**: 45 diagrams needed. Priority order? Are some more important than others?

### Intentionally Excluded Content

- 4 IPhO note photos out of scope: cover page, lab safety, exam preparation tips, notebook boundary marker
- Lab sections and error propagation (per foreword: "no lab sections")
- Test-taking strategies (per foreword)

---

## Pre-Casey Sync (DONE)

- [x] Physics correctness pass across all IPhO-note chapters
- [x] Consistency pass (notation, units, formatting) across all chapters
- [x] Status summary (above)

## Phase 3: Polish & Publish (after Casey sync)

### Diagrams
- [ ] Replace ~45 remaining diagram placeholders with TikZ
- [ ] Review all AI-generated diagrams for correctness (remove "Needs Review" callouts after sign-off)

### Practice Problems
- [ ] Add 7-Ds practice problems to: thermo, electrostatics, circuits, optics, relativity

### Content
- [ ] Write co-author preface (Jonathan + Casey)
- [ ] Unify tone across chapters (pending Casey feedback)
- [ ] Add cross-references between chapters where natural

### Infrastructure
- [x] ~~Add PDF rendering to CI~~ (done — CI renders both HTML and PDF, sidebar has download link)

---

## Open Questions
1. Any topics Casey wants to add beyond what's in the IPhO notes?
2. Error propagation — considered under lab stuff that we're not covering?
3. Should the math chapter get more physics-motivated examples?
