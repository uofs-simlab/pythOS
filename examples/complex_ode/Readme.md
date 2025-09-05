# Complex ODE problem


## Problem Description

```math
\begin{aligned}
\frac{\rm{d}u}{\rm{d}t} = iu + 0.1u - 0.1u^3,
\end{aligned}
```
which may be solved as a system of real ODEs
```math
\begin{aligned}
\frac{\rm{d}x}{\rm{d}t} & = - y+  0.1x  +0.3xy^2 -0.1x^3,  \\
\frac{\rm{d}y}{\rm{d}t} & = \hspace*{1.5ex}x+0.1y -0.3x^2y + 0.1y^3.
\end{aligned}
```
The problem is solved with initial condition $`u(0) = 0.1`$, and $` t \in [0,100] `$


## Contents

The complex_ode.py file demonstrates solving the problem as both a real and complex system using the fractional step solvers