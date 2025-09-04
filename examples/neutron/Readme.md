# Neutron Transport Equation

This is one of the examples included in M. M. Moayeri, V. Guenter, S. Wei, and R. J. Spiteri, "pythOS: An Operator-Splitting Library in Python"

## Problem Description

```math
\begin{aligned}
\frac{1}{v}\frac{\partial\phi}{\partial t} & = \nabla\cdot(D\nabla\phi) - \Sigma_a\phi + S(t,\mathbf{x})
\end{aligned}
```
with boundary conditions
```math
\begin{aligned}
  \phi &= \bigl.\phi_D(t,\mathbf{x}) = 0 \bigr|_{\mathbf{x}\in\partial\Omega_D}, \\
  \Bigl.(D\nabla\phi)\cdot\mathbf{n}\Bigr|_{\mathbf{x}\in\partial\Omega_N} &= \Bigl.\phi_N(t,\mathbf{x}) = 0 \Bigr|_{\mathbf{x}\in\partial\Omega_N}, \\
\partial\Omega_D &= \left\{(0,y),(5,y)\right\}_{y\in[0,5]} \\
\partial\Omega_N &= \left\{(x,0),(x,5)\right\}_{x\in[0,5]}
\end{aligned}
```

The spatial domain is defined as $`\mathbf{x} \in [0,5]^2 `$, and the problem is solved over $` t  \in [0,10] `$

Initial conditions are defined as 
```math
\begin{aligned}
\phi(\mathbf{x}) &= 0
\end{aligned}
```

The source term is defined as 
```math
S(t,\mathbf{x}) = 10\left(\frac{\omega}{v}\left(bx-x^2\right)\cos\left(\omega t\right)+\Sigma_a\left(bx-x^2+2D\right)\sin\left(\omega t\right)\right),
```
with constants defined such that the problem has an exact solution
```math
\phi(t,x) = A\sin(\omega t)(bx-x^2).
```
The domain is discretized with continuous Galerkin FEM with quadratic elements on a 2D rectangular mesh implemented through Firedrake.

## Contents

The neutron.py file demonstrates the solution using fractional step methods (specifically shown is Strang) across multiple $` \Delta t`$ to reproduce the convergence results from the paper