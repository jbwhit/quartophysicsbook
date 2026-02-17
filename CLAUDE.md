# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

"Think Like A Physicist" — a terse undergraduate physics textbook by Casey Handmer and Jonathan Whitmore. Quarto book that renders to both a website and PDF. Source material includes Casey's "Think Like A Physicist" PDF and 169 handwritten IPhO (International Physics Olympiad) training note photos in `images/IPhO notes/`.

Scope: physics thinking and problem-solving only. No lab sections, no test-taking content.

## Workflow

Commit early and often, and push to GitHub after each logical chunk of work.

### TODO Hygiene

`TODO.md` tracks current and upcoming work. When a logical section (e.g., an entire Phase) is fully complete, move its items from `TODO.md` into `DONE.md` to keep the TODO file clean and focused on what's left.

## Build & Development Commands

```bash
# Preview site locally (live reload)
quarto preview

# Render full site (output goes to _site/)
quarto render
```

## Architecture

Quarto book project at repo root. Chapter files are `.qmd` at top level, organized into parts in `_quarto.yml`. Old attempts archived in `_archive/` (gitignored). Source images in `images/`. Output goes to `_book/` (gitignored).

## Content Conventions

### The 7 Ds Problem Template

Every worked problem MUST follow the 7 Ds (and the little s) structure. Each step must appear as a subsection, even if the content is brief. If a step is not applicable, write a short explanation (e.g., "Not needed for this problem" or "Diagram not needed — purely algebraic").

```markdown
## Problem Title {.unnumbered}

Problem statement goes here.

::: {.callout-warning}
## TODO: Diagram
Description of what the diagram should show. Include key elements, labels, and coordinate system.
:::

### Diagram {.unnumbered}
[Diagram or placeholder callout. Should be big — 2/3 of a page. Load the problem into your GPU.]

### Directions {.unnumbered}
[Define positive/negative directions. State coordinate system choice.]

### Definitions & Given Data {.unnumbered}
[Table or list of ALL variables with descriptions and units. Put everything on the page.]

### Diagnosis {.unnumbered}
[What type of problem is this? Identify the relevant physics: conservation laws, force laws, etc.]

### Derivation {.unnumbered}
[Write fundamental equations. Express diagnosis in symbols. Need as many equations as unknowns. Check dimensions as you go.]

### Determination (D'algebra) {.unnumbered}
[Algebraic manipulation to isolate the answer. Box the final result with $$\boxed{...}$$]

### Dimensions {.unnumbered}
[Verify dimensions of the answer. Check limiting cases / sanity check. Use smiley or checkmark for each case that makes sense.]

### Substitution {.unnumbered}
[Only if needed. Rough numerical calculation by hand. Check units. Include error estimate.]
```

### Diagram Placeholders

When a diagram is needed but not yet created, use a warning callout describing what the diagram should contain:

```markdown
::: {.callout-warning}
## TODO: Diagram
Description of diagram contents, key elements, labels, etc.
:::
```

### IPhO Note Photos

The handwritten IPhO note photos in `images/IPhO notes/` are reference material only — do NOT embed them in the book. Transcribe their content into proper Quarto markdown with LaTeX math.


## Deployment

Push to `main` triggers `.github/workflows/publish.yml`:
1. Quarto renders the site to `_site/`
2. `actions/deploy-pages` deploys to GitHub Pages

GitHub Pages source must be set to "GitHub Actions" (not "Deploy from a branch").

## Theming

Dual light/dark mode using flatly (light) and darkly (dark) Bootstrap themes, with custom SCSS overrides.
