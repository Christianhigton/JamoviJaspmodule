# Jamovi JASP Enhanced Modules

This workspace contains jamovi module development for bringing selected JASP-style analyses into jamovi with a native jamovi interface.

Current active project:

- `jaspEnhancedAnova`: focused build exposing only:
  - JASP Enhanced ANOVA
  - JASP Enhanced Repeated Measures ANOVA

Generated module artifact:

- `jaspEnhancedAnova_0.1.0.jmo`

The regression scaffold is present but not part of the current active deployment.

## Project Scope

The current goal is to stabilize the two implemented ANOVA analyses before re-enabling additional analysis menus. Disabled ANOVA-family analysis shells are kept under `jaspEnhancedAnova/disabled/` for later work.

## Build

See [DEPLOYMENT.md](DEPLOYMENT.md) for exact build and install commands.

## Deployment Page

A static project page is available at [docs/index.html](docs/index.html). It is suitable for GitHub Pages if this repository is published and Pages is configured to serve from the `docs/` directory.

## License and Attribution

This project follows GPL (>= 2), matching the upstream JASP ANOVA module license declared in its DESCRIPTION file.

Upstream reference:

- https://github.com/jasp-stats/jaspAnova

Please cite both jamovi and JASP when using adapted functionality.
