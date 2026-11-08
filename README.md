# CylinderFlow: Interactive 2D CFD Simulator Using the Lattice Boltzmann Method

A self-contained, browser-based simulator of two-dimensional viscous flow past a circular cylinder. It uses a **D2Q9 lattice Boltzmann solver (BGK collision)** written in vanilla JavaScript. The app is one file, `index.html`, with nothing to install, no libraries, no server and no internet connection needed.

## Run it

1. Download `index.html` (or clone the repo).
2. Double-click it to open it in Chrome, Edge or Firefox.
3. Click a preset, e.g. **Intermediate · 100**.

## What it does

- **Numerical solver:** D2Q9 lattice, BGK collision, fused pull-streaming, double precision.
  - Boundaries: half-way bounce-back cylinder (no-slip), equilibrium velocity inlet, zero-gradient outlet with an optional absorbing layer, free-slip side walls.
  - Stability: runs are checked continuously and stopped with an explanation if they diverge.
- **Forces:** drag and lift from the momentum-exchange method. The Strouhal number comes from lift zero crossings, and the recirculation length is measured on the centreline.
- **Visualisation:** velocity magnitude, vorticity, pressure coefficient and instantaneous streamlines. Tracer particles, a hover probe and a snapshot filmstrip show the wake over time.
- **Analysis:** live metrics, plus charts of Cd(t*), Cl(t*), mean velocity, Cd vs Re and surface Cp(θ). You can save cases and compare them.
- **Export:** time series as CSV, full results as JSON, the instantaneous field as CSV, and images as PNG.
- **Built-in documentation:** Theory, About & Methodology (abstract, governing equations, LBM formulation, boundary conditions, stability, results, limitations, references) and an in-app **Verification** suite.

## Verification

The solver is checked against analytical results. The suite runs in the app (Verification tab) and from the command line:

```
node scripts/run_tests.js
```

**Result: 24 / 24 checks pass.** They include:

| Check | Result |
|---|---|
| D2Q9 weights and moment isotropy, equilibrium moments | exact to round-off |
| Mass and momentum conservation (collision + streaming) | drift ~1e-15 |
| Viscosity from shear-wave decay vs ν = (τ − ½)/3, τ = 0.6/0.8/1.2 | error < 0.1 % |
| Poiseuille flow between bounce-back walls | L2 error 0.07 % |
| Momentum-exchange force = fluid momentum lost each step | ~1e-12 |
| Uniform flow preserved by inlet/outlet/free-slip boundaries | ~1e-14 |
| Instability detection | stops safely |

The finished `index.html` was also exercised offline in headless Chrome: 32 / 32 interface checks passed (controls, presets, reset, case comparison, all exports, full screen, persistence, responsive layout, no console errors). See `results/ui_test_results.json`.

**Screenshots:** open `index.html` and press a preset; the vortex street at Re = 100 forms within about 40 convective time units.

## Results

All results below were produced by the same solver code, run in Node.js (`scripts/reference_run.js`, raw output in `results/reference_runs.jsonl`). Every case uses U = 0.08 (lattice units), 12.5 % blockage and the cylinder 4D from the inlet.

| Re | Grid | D (cells) | Mean Cd | Cl amplitude | St | Lr/D | Wake |
|---|---|---|---|---|---|---|---|
| 20 | 400×160 | 20 | 2.659 | 0.0004 | — | 1.011 | steady |
| 40 | 400×160 | 20 | 1.984 | 0.004 | — | 2.328 | steady |
| 100 | 400×160 | 20 | 1.671 | 0.427 | 0.1825 | — | shedding |
| 180 | 400×160 | 20 | 1.660 | 0.921 | 0.2056 | — | shedding |
| 100 | 300×120 | 15 | 1.686 | 0.427 | 0.1805 | — | shedding |
| 100 | 600×240 | 30 | 1.618 | 0.446 | 0.1841 | — | shedding |

**How to read these numbers**

- **The flow regimes come out right:** a steady separated wake at Re 20 and 40, and periodic vortex shedding at Re 100 and 180.
- **Recirculation lengths are close to the published unconfined values** (Dennis & Chang 1970: about 0.94 at Re 20 and 2.35 at Re 40).
- **Drag and Strouhal number come out higher than unconfined benchmarks** (e.g. Park, Kwon & Choi 1998 at Re 100: Cd ≈ 1.33, St ≈ 0.165). The main cause is domain confinement. On a larger domain (800×480, 4.2 % blockage, 10D upstream) the Re 20 drag fell from 2.66 to 2.20.
- **Grid refinement from D = 15 to D = 30 changes Cd by about 4 %.**

No experimental validation has been performed. The published values are quoted only as context.

## Limitations

The model is 2-D and laminar, and represents the cylinder as a staircase of grid cells. Other limitations:

- Strong confinement (12.5 % blockage) and simple inlet/outlet boundary conditions.
- Weakly compressible (Mach ≈ 0.14–0.23 locally).
- BGK collision limits the reachable Re for a given resolution.

The full list, and the planned improvements, are in the app's About & Methodology panel.

## Repository layout

```
index.html                 the complete application (solver + UI + docs, single file)
src/solver.js              solver source (identical to the code inlined in index.html)
src/tests.js               verification suite (identical to the in-app suite)
scripts/run_tests.js       run the verification suite with Node.js
scripts/reference_run.js   reproduce a row of the results table
results/                   raw reference-run output and UI test log
```

## References

- Bhatnagar, Gross & Krook (1954), *Phys. Rev.* 94, 511–525.
- Qian, d'Humières & Lallemand (1992), *Europhys. Lett.* 17, 479–484.
- Chen & Doolen (1998), *Annu. Rev. Fluid Mech.* 30, 329–364.
- Krüger et al. (2017), *The Lattice Boltzmann Method: Principles and Practice*, Springer.
- Ladd (1994), *J. Fluid Mech.* 271, 285–309.
- Mei, Yu, Shyy & Luo (2002), *Phys. Rev. E* 65, 041203.
- Dennis & Chang (1970), *J. Fluid Mech.* 42, 471–489.
- Park, Kwon & Choi (1998), *KSME Int. J.* 12, 1200–1205.
- Williamson (1996), *Annu. Rev. Fluid Mech.* 28, 477–539.

---
Author: Pala Geya Sri Sai Harshit, B.Tech Engineering and Computational Mechanics, IIT Delhi.
