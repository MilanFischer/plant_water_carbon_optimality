# Nakad 2023 AFM Explorer (Shiny)

Interactive Shiny app implementing a **simplified** version of the stomatal optimality framework in **Nakad et al. (2023)** to explore:
- optimal stomatal conductance response to VPD,
- photosynthesis and transpiration responses under an optimality constraint,
- **AFM** (atmospheric feedback mechanism) onset diagnostics,
- optional **marginal WUE** (mWUE = dfc/dfe) when λ is inferred from a reference ci/ca.

> Note: This is an educational / exploratory tool. It reproduces core algebraic relationships used in the paper’s simplified analytical development, not a full land-surface or full biochemical (Farquhar) model.

## Features

- **Two λ modes**
  - **Diagnostic mode**: infer λ by inverting Eq. 20 at a chosen `VPDref` and `ci/ca(ref)`
  - **Forward mode**: set λ directly
- Plots vs VPD:
  - `gCO2`, `gH2O = a·gCO2`
  - `fc` (photosynthesis), `fe` (transpiration proxy)
  - `ci/ca` (Eq. 20 in open-stomata regime; `ci/ca → 1` when `g=0`)
  - `WUE = fc/fe` and closed-form Eq. 21
  - **mWUE = d(fc)/d(fe)** (only when diagnostic λ mode is enabled), with λ shown as a reference line
- **AFM diagnostics**
  - theory onset `Dafm` (shown as vertical dashed line)
  - numeric `fe` peak and derivative zero-crossing (shown as dotted line)

## Model summary (simplified)

The app uses a compact set of relationships consistent with the simplified analytical framework:
- Fick-type and transpiration proxy: `fe ≈ a · gCO2 · D`
- A demand form for photosynthesis: `fc(g)`
- Optimal stomatal conductance: `g*(D)`
- Optimal `ci/ca(D)` and WUE expressions
- AFM onset estimate `Dafm`

Where:
- `D = VPD / Patm` is used (dimensionless mol/mol) for unit consistency
- `a ≈ 1.6` is the diffusivity ratio relating H2O and CO2 conductances

## Getting started

### Requirements
- R (>= 4.0 recommended)
- Packages: `shiny`, `ggplot2`

### Install packages
```r
install.packages(c("shiny", "ggplot2"))
