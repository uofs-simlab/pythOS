# Additive Runge Kutta integrators

The additive Runge Kutta integrator may be used with any number of operators. 

The relevant integrator function is `ark_solve`. It takes some basic arguments to define the problem and method, and optional arguments may be used to control output and integrator behaviour

## Basic Usage

The `ark_solve` function takes as inputs:

- the operators to use
     	 
	A list of operators (see [Operators](../README.md#operators)) in the order you wish to use them.

- dt

	the time step to advance by, which may be an integer or a float

- y0

	the initial value(s) to use for y (see [Initial Condition](../README.md#initial-condition))

- t0

	the time value to start at
	
	If using the finite element version, this must be of type Constant
- tf

	the time value to advance to

- methods
	a list of `Tableau` or `EmbeddedTableau`, one for each operator.

	See [Method Definition](#method-definition) for further details on the definition of the method.
	
	Note any extra tableaus are ignored.  If there are fewer tableau than operators, no computation is performed.

## Method Definition

The method is defined by a list of `Tableau` or `EmbeddedTableau`.

A `Tableau` is a class provided in the file `butcher_tableau.py`, which takes as inputs:
- c : a 1 dimensional numpy array for c
- a : a 2 dimensional numpy array for a
- b : a 1 dimensional numpy array for b


The `EmbeddedTableau` class is provided in the `butcher_tableau.py` file, which takes as inputs:
- c : a 1 dimensional numpy array for c
- a : a 2 dimensional numpy array for a
- b : a 1 dimensional numpy array for b
- b\_aux : a 1 dimensional numpy array with the embedded method
- order : the lesser of the orders of the two methods, used in step size calculations

Depending on the structure of the methods, pythOS will use an explicit, diagonally implicit, or fully implicit additive Runge Kutta integrator.

If all supplied tableaus are of type `EmbeddedTableau`, the adaptive solver options is enabled. If any tableaus are of type `Tableau`, any adaptive information from `EmbeddedTableau` or tolerance specifications is ignored.

Note that the integrator assumes that all methods are the same size (i.e. have the same number of underlying stages)

## Optional Arguments

### File Output

- fname: a filename to save intermediate results to.  
  
	 When using the finite element capabilities of firedrake this is a 
	 .h5 file accessible through the `CheckpointFile` interface from 
	 `firedrake`. The file has attributes `/times/idx/` that store the time 
	 of each entry and the attribute `/time/last_idx` indicates the last 
	 valid index. The `Function`s are stored using corresponding indices. 
	 Otherwise this is a .csv file containing time in the first entry of each 
	 line, and the solution vector in the remaining entries

- save\_steps: the number of intermediate steps to save.  
  
	 The default is to save steps at delta\_t time interval if a filename
	 is provided.
	 Note that if an adaptive method is used, the saved steps may not be evenly spaced

### Adaptive Methods

- rtol: The relative tolerance to use if all supplied tableau are `EmbeddedTableau`
- atol: The abolute tolerance to use if all supplied tableau are `EmbeddedTableau`

### Dynamic Linearization

- jacobian: A function to compute a jacobian when using dynamic linearization

	If this option is used, the provided operator functions must take as arguments (t, y, J)

When this is used, the jacobian function is called at the start of each step, and held constant for each substep within.  This jacobian is passed to each operator function, which may be uses to compute the values of each operator

### Boundary Condition

- bc: Any boundary conditions to apply when using firedrake

### Controlling Solver Behaviour

An optional dictionary of additional arguments may be supplied through the `solver_parameters` argument.

If it is supplied, it is used as optional arguments for `scipy.optimize.root` if the operators are defined as functions, or as optional arguments for firedrake's `solve` function if operators are defined as `Form`s. 

For the firedrake finite element integrator, the `solve` function is applied to determine the new stage value (referred to as y), and then the stage derivatives (k=f(y)). The optional arguments are only applied to determining the stage values y.
