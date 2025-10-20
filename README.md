# formin_kinetic_model

Contains code to compute polymerization rates for binding sites (PRMs) on FH1 domains and replicate figure 4 of Bogue et al., *Formin N-terminal dimerization can impact F-actin assembly via polymer mechanics of the FH1 domain*, publication pending.

## Contents

* ``src/``: MATLAB scripts and functions to generate a simple example figure.
* ``runs/``: where the output figure is saved; includes an example file.

## Quick Start

#### Set-up

1. Clone this repository to your machine
2. Install MATLAB (The MathWorks, Inc.) if not already installed

#### Generating Example plot

1. Either open the MATLAB GUI or run MATLAB via command-line.
2. Navigate to the ``src/`` directory.
3. Run `makeKineticHeatmaps.m`

This will generate plots of of polymerization rate ratios (dimer/nondimer) for sweeps across possible probability density and occlusion probability ratios for 3 overall parameter regimes (capture rate-limiting, capture and delivery competing, and delivery rate-limiting). The plots will be saved as kineticHeatmaps.png in the ``runs/`` directory.

The output should be the following plot:
![kinetic sweep heatmaps for 3 parameter regimes](runs/kineticHeatmaps_example.png)
