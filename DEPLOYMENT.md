# Deployment

These commands build and install the current `jaspEnhancedAnova` development module.

## Active Module

Only these analyses are currently enabled:

- JASP Enhanced ANOVA
- JASP Enhanced Repeated Measures ANOVA

Disabled analyses are retained in `jaspEnhancedAnova/disabled/`.

## Prepare

```sh
cd /Users/christianhigton/JamoviJaspmodule/jaspEnhancedAnova
'/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/node/node-darwin/bin/node' \
  '/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/jmvtools/node_modules/jamovi-compiler/index.js' \
  --prepare '.' \
  --home '/Applications/jamovi.app' \
  --assume-app-version '2.7.0'
```

## Build `.jmo`

```sh
cd /Users/christianhigton/JamoviJaspmodule/jaspEnhancedAnova
'/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/node/node-darwin/bin/node' \
  '/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/jmvtools/node_modules/jamovi-compiler/index.js' \
  --build '.' \
  --home '/Applications/jamovi.app' \
  --assume-app-version '2.7.0' \
  --jmo '/Users/christianhigton/JamoviJaspmodule/jaspEnhancedAnova_0.1.0.jmo'
```

## Manual Local Install

The jamovi binary on this machine does not accept `--install`, so install by unpacking into the user module directory:

```sh
cd /Users/christianhigton/JamoviJaspmodule
rm -rf /private/tmp/jaspEnhancedAnova-install
unzip -q jaspEnhancedAnova_0.1.0.jmo -d /private/tmp/jaspEnhancedAnova-install
ditto /private/tmp/jaspEnhancedAnova-install/jaspEnhancedAnova \
  "/Users/christianhigton/Library/Application Support/jamovi/modules/jaspEnhancedAnova"
```

Fully quit and reopen jamovi after installing.

## Verify Installed R Package

```sh
Rscript -e 'library(jaspEnhancedAnova, lib.loc="/Users/christianhigton/Library/Application Support/jamovi/modules/jaspEnhancedAnova/R"); cat("loaded\n")'
```

## Re-enable Disabled Analyses Later

To re-enable an analysis:

1. Move its `.a.yaml`, `.u.yaml`, and `.r.yaml` files from `jaspEnhancedAnova/disabled/jamovi/` back into `jaspEnhancedAnova/jamovi/`.
2. Move its `.b.R` file from `jaspEnhancedAnova/disabled/R/` back into `jaspEnhancedAnova/R/`.
3. Run prepare and build again.
4. Reinstall the rebuilt `.jmo`.
