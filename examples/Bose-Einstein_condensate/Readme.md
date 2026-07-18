# 2-dimensional coupled Gross-Pitaevskii equation for Bose-Einstein condensate 

A 2D 2-component Bose-Einstein condensate with Rabi coupling is simulated using the coupled Gross-Pitaevskii (non-linear Schrödinger) equations. 

## Problem Description
$$
\begin{aligned}
i\frac{\partial}{\partial t}\Psi_1 &= \frac{1}{2}\bigg[-\nabla^2 + y^2 + R_0^2|\Psi_1|^2 + (1+\gamma)R_0^2|\Psi_2|^2 \bigg] \Psi_1 - \omega\Psi_2, \\
i\frac{\partial}{\partial t}\Psi_2 &= \frac{1}{2}\bigg[-\nabla^2 + y^2 + (1+\gamma)R_0^2|\Psi_1|^2 + R_0^2|\Psi_2|^2 \bigg] \Psi_2 - \omega\Psi_1.
\end{aligned}
$$

The domain is a long ribbon, $x\in[-1.72,1.72]$, $y\in[-40.96,40.96]$. The problem is solved for $t\in[0,1]$. We use $R_0 = 30$ which controls the non-linear component couplings, $\gamma = 0.01$ which is a repulsion parameter, and $\omega = 1.0$ which sets the Rabi frequency. The condensates are initially phase-separated with a small sinusoidal perturbation given to their shared interface.

Spatial step size used is $dx = dy = 0.04$, and timestep size $dt = 0.0004$.

Operator splitting is used along with analytical solutions for the different operators. A pseudospectral method is used for the Laplacian.

## Contents

This code uses the `pythOS` library.

The file `BEC.py` creates .csv files containing the solution at time intervals of $\Delta t = 0.1$.

The file `BEC_timing_files.py` creates .csv files containing the solution at time intervals of $\Delta t = 0.1$, and a .txt file recording CPU runtime. To instead test multiple different $dt$ and record final runtimes and errors (calculated from a reference solution, which requires generation), set the variable 'deltat' equal to the variable 't_end'.

The file `plot_2D.py` creates 2D plots of the component densities and the relative phase of the condensate, using data generated from `BEC.py`.

The file `plot_error_CPU.py` creates a plot of data generated from running tests of multiple different $dt$ using `BEC_timing_files.py`.
