# Stiff Brusselator problem

This folder includes files demonstrating the solution of the stiff brusselator problem using GARK, multirate and multirate_infinitesimal solvers

## Problem Definition

```math
\begin{aligned}
u_t &= \alpha_u u_{xx} + \rho_u u_x + a - (w + 1)u + u^2v \\
v_t &= \alpha_v v_{xx} + \rho_v v_x + wu - u^2v \\
w_t &= \alpha_w w_{xx} + \rho_w w_x + (b-w)/\epsilon - wu
\end{aligned}
```

with $`\alpha_u = \alpha_v = \alpha_w = 10^{-2}`$, $`\rho_u=\rho_v=\rho_w = 10^{-3}`$, $`a=0.6`$, $`b=2`$ and $`\epsilon = 10^{-3}`$

The domain is discretized into 201 uniformly spaced points in the domain $`x \in [0, 1]`$

## Solvers demonstrated - details of methods are contained within the files

- `stiff_brusselator_gark.py`: GARK solver with second order method
- `stiff_brusselator_multirate.py`: Multirate GARK solver with second order method
- `stiff_brusselator_multirate_infinitesimal.py`: Multirate Infinitesimal methods with both IRK and IMEX slow stages
