# Stepsize selection for time-varying online feedback optimization using spectral information

## Authors

Fritz A. Engeln\*, Marta Zagórowska\*, Lukas Ortmann\*\* , Jan-Willem van Wingerden\*

\* Delft Center for Systems and Control, Delft University of Technology, Delft, The Netherlands

\*\* Institute for Communication Systems, Eastern Switzerland University of Applied Sciences, Rapperswil, Switzerland

Corresponding author: Fritz A. Engeln
Email: f.a.engeln@tudelft.nl
## Abstract

Time-varying convex optimization problems commonly arise in feedback optimization controllers, where changes in system dynamics, operating conditions, or external disturbances can cause the minimizer to evolve over time. We analyze the tracking performance of a scaled gradient algorithm using contraction theory and derive a stepsize that minimizes an upper bound on the tracking error. We then extend the analysis to exploit knowledge of the Fourier spectrum of the time-varying minimizer to determine a stepsize that further reduces the tracking error. The proposed results are demonstrated on a wind farm feedback optimization problem, where the frequency of the minimizer is inferred from changes in wind direction. Numerical results show that incorporating spectral information into stepsize selection improves tracking performance compared with a stepsize designed without prior spectral information.

## Overview

This repository contains the code accompanying the paper:

**"Stepsize selection for time-varying online feedback optimization using spectral information"**

The repository provides all scripts required to reproduce the numerical example presented in the paper.
<p align="center">
  <img src="WindTurbines/figures/turbinesSetup.png" alt="Two-tank system" width="800">
</p>

## Reproducing the numerical results

The numerical results presented in this work can be reproduced using the code provided in this repository. The results are divided into two examples:

1. **Illustrative example with a time-varying reference**
2. **Wind farm control example**

### 1. Illustrative example: Time-varying reference

All files required to reproduce the illustrative example are located in the `IllustrativeExample/` folder.

The following scripts can be used to reproduce the results:

- **Run the OFO controller:** Execute `IllustrativeExample/OFOController/main` to run the online feedback optimization (OFO) controller with different step sizes.
- **Calculate the optimal step size and the function \(J\):** Execute `IllustrativeExample/CalcStepsize` to calculate the optimal step size and evaluate the function \(J\) based on spectral information.
- **Plot the results:** Execute `IllustrativeExample/plotErrors` to generate plots of the resulting tracking errors.

### 2. Wind farm control example

All files required to reproduce the wind farm control example are located in the `WindTurbines/` folder.

The following scripts can be used to reproduce the results:

- **Run the OFO controller:** Execute `WindTurbines/mainOFO` to run the OFO controller with different step sizes.
- **Calculate the optimal step size and the function \(J\):** Execute `WindTurbines/calcStepsize` to calculate the optimal step size and evaluate the function \(J\) based on spectral information of the wind direction.
- **Plot the tracking error:** Execute `WindTurbines/plotErrors` to visualize the resulting tracking error.
- **Plot the wind direction:** Execute `WindTurbines/plotWind` to visualize the wind direction used in the wind farm control example.

<!--
## Citation

If you use this code in your research, please cite:

F. A. Engeln, M. Zagórowska, and J.-W. van Wingerden, *Time-varying feedback optimization: Convergence bounds and stepsize selection*
--!>
