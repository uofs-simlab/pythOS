# 1-dimensional monodomain problem

This is one of the examples included in M. M. Moayeri, V. Guenter, S. Wei, and R. J. Spiteri, "pythOS: An Operator-Splitting Library in Python"

## Problem Description

```math
\begin{aligned}
    \chi C_m\frac{\partial v}{\partial t}   &= \frac{\lambda}{1+\lambda}\nabla\cdot(\sigma_i\nabla v) - \chi \mathbf{I}_{\rm{ion}}(s,v),\\
    \frac{\partial s}{\partial t}  &= \mathbf{g}(s,v)
\end{aligned}
```
with homogeneous Neumann boundary conditions.

In this example, the FHN model is used:
```math
\begin{aligned}
\mathbf{I}_{\rm{ion}}(s,v)  &= -k\left(v- v_{r} \right) \left( s+\left( v-v_{th}\right)\left(v-v_{p}\right)\right),\\
    g(s,v) &=L\left(v-v_{r}\right)-bs,
\end{aligned}
```

The spatial domain is defined as $`x \in [0,1] `$, and the problem is solved over $` t  \in [0,15] `$

Initial conditions are defined as 
```math
\begin{aligned}
v_0(x) &= \begin{cases}
      -85+100\cos^2(\pi x/0.2), & x\in [0,0.1], \\
      -85, & \mathrm{otherwise},
       \end{cases} \\
s_0(x) &= 0
\end{aligned}

```

The domain is discretized with finite differences with a uniform spatial grid containing 131 points

## Contents

The 1D_FHN_monodomain.py file demonstrates the solution using fractional step methods across multiple $` \Delta t`$