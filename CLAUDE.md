# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview


## Workflow

Commit early and often, and push to GitHub after each logical chunk of work.

## Build & Development Commands

```bash
# Preview site locally (live reload)
quarto preview

# Render full site (output goes to _site/)
quarto render
```

## Architecture

TODO

## Content Conventions


## Deployment

Push to `main` triggers `.github/workflows/publish.yml`:
1. Quarto renders the site to `_site/`
2. `actions/deploy-pages` deploys to GitHub Pages

GitHub Pages source must be set to "GitHub Actions" (not "Deploy from a branch").

## Theming

Dual light/dark mode using flatly (light) and darkly (dark) Bootstrap themes, with custom SCSS overrides.
