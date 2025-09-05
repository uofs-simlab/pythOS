# 2-dimensional advection-diffusion-reaction problem


## Problem Description

```math

\begin{aligned}
u_t = - \alpha \left(\nabla \cdot u\right) + \epsilon
  \left(\nabla^2 u\right) + \gamma u \left(u-1/2\right)\left(1-u\right)
\end{aligned}
```
with homogeneous Neumann boundary conditions.


The spatial domain is defined as $`x \in [0,1]^2 `$, and the problem is solved over $` t  \in [0,0.1] `$, and with parameters $\alpha = -10$, $\epsilon=1/100$, $\gamma = 100$

Initial conditions are defined as 
```math
\begin{aligned}
  u\left(x, y, 0\right) = 256\left(xy\left(1-x\right)\left(1-y\right)\right)^2 + 0.3
\end{aligned}

```

The domain is discretized with a uniform spatial grid with $` \Delta x = \Delta y = 1/40 `$, using the continuous Galerkin finite element method with
quadratic elements as implemented in Firedrake.

## Contents

The 2D_ADR.py file demonstrates an order convergence study for three and four-splittings using fractional step methods.