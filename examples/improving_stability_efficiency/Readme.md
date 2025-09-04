# Improving the stability and efficiency of high-order operator-splitting methods

This repository contains the source codes for the experiments presented in the paper S. Wei, V. Guenter, and R.J. Spiteri, "Improving the stability and efficiency of high-order operator-splitting methods"

## Problem Definition

```math
\begin{aligned}
\chi C_m v_t + \chi \mathbf{I}_{\rm{ion}}(\mathbf{s},v) &= \nabla \cdot (\sigma \nabla v), \\
\mathbf{s}_t &= \mathbf{f}(t,\mathbf{s},v), 
\end{aligned}
```

The domain is $`x \in [0, 2], y \in [0, 0.7], z \in [0, 0.3]`$, solved over a time domain of $`t \in [0, 40]`$

The domain is discretized with a uniform grid with $\Delta x = \Delta y = \Delta z = 0.05$

The spatial derivatives are approximated with second-order centered differences.

## Contents
The file `Niederer.py` contains the code implementing the problem.  To run an experiment, change the dt, alpha and experiment variables as required.

The code uses the `pythOS` library.  For details on how to set up this library, see the main [README](https://github.com/uofs-simlab/pythOS/blob/main/README.md)

The file `Niederer_cell_ml.py` is a later addition that solves the Niederer problem using the CellML definition of the ten Tusscher & Panfilov model. This demonstrates an easy way to change the cell model for other cardiac experimentation.