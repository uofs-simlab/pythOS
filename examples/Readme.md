# Examples for using the pythOS library

This directory contains examples of setting up a problem, and using a variety of solver options within the library.

Included in this directory are:

- 1D_monodomain: A 1D monodomain problem using FHN cell model, solved using the fractional step solver
- advection_diffusion_reaction: A 2D advection-diffusion-reaction problem, discretized using Finite element method through Firedrake, and demonstration of an order convergence study for 3 and 4 splittings with fractional step methods
- complex_ode: A complex ODE demonstrating solution as a real system or complex ODE with the fraction step solvers
- improving_stability_efficiency: Niederer benchmark problem using the fractional step solver (from the paper S. Wei, V. Guenter, and R.J. Spiteri, "Improving the stability and efficiency of high-order operator-splitting methods"). It also includes a file to demonstrate integration of CellML models for other cardiac problems.
- neutron: A neutron transport equation, discretized using Finite element method through Firedrake, with an exact solution and demonstration of order calculations for fractional step methods.
- stiff_brusselator: The stiff brusselator problem, solved using GARK, multirate and multirate_infinitesmial solvers
- tests: a series of simple test runs for the various solvers that use the brusselator ODE problem, which may be useful to ensure that the library is correctly set up for a new user.
