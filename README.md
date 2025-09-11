# pythOS
	
## Setup

Requirements:

- numpy
- scipy

Optional:

- firedrake - for finite element capabilities
- Irksome - for additional finite element capabilities
- SUNDIALS (v7) - to use any of the sundials integrators for a subintegration

## Usage

The relevent code must be on the `PYTHONPATH` prior to import.

On Unix systems, this can be achieved by running

```
export PYTHONPATH=$PWD:$PYTHONPATH

```
from the root directory of pythOS.

On Windows systems, this can be acheived by running 

```
set PYTHONPATH=%CD%;%PYTHONPATH%
```

### Using SUNDIALS:

If using SUNDIALS, prior to use, the wrapper library must be built.

On Unix systems, a Makefile is provided to build the wrapper.  Both the wrapper and the SUNDIALS lib folder must be in the `LD_LIBRARY_PATH`.

On Windows systems, the build has been tested using the Visual Studio build tools. The pythOS wrapper code expects the SUNDIALS library to be installed in a dependencies/install_sundials sub-directory.  Following this, the wrapper can be built with

```
cl /LD sundials_wrapper.c -I dependencies\install_sundials\include
```

Then the sundials wrapper can be used in pythOS without any additional setup.

### Basic structure of arguments

All the integrators require some basic information:

#### Operators

Depending on the problem being solved and the integrators in use, the operators may take a few different formats:
- functions with arguments (time, y vector) that returns another vector containing dydt
- functions with arguments (time, y vector, Jacobian) that returns another vector containing dydt when the dynamic linearization option is enabled.
- functions with arguments (time, y `Function`) that returns a `Function` representing dydt (for use with some options with the fractional step integrator only)
- `Forms` representing the value of the operator (dydt) when using Firedrake
- (`Form`, boundary condition) when a boundary condition is required
- Irksome `Forms` which must contain any required `Dt` terms when using an Irksome integrator.
- functions with arguments (time, delta t, y vector) that returns a vector for an analytical solution
- functions with arguments (time, delta t, y `Function`) that returns a `Function` for an analytical solution
- functions with arguments (time, y vector, dydt vector) that returns a residual vector for a solution using IDA

Some of these options only work with some of the integrators and sub-integrator choices.  This will be further detailed in the sections for the relevant section

#### Initial Condition

This has two structures:
- a numpy array (which should be 1-dimensional)
- a Firedrake `Function` is required when finite element capabilities are used through Firedrake

#### Time information

Each integrator requires an initial time, final time and time step (delta t). These are all numeric values.  They may also be of type `Constant` if Firedrake is in use.
Depending on the integrators used, the initial time and delta t value may need to be of type `Constant` for Firedrake.

### Integrator Options

To jump to the basic documentation for each integrator, use the links here:
- [Fractional Step integrator](#fractional-step-integrator)
- [Additive Runge Kutta integrator](#additive-runge-kutta-methods)
- [General structure additive runge kutta methods](#general-structure-additive-runge-kutta-methods)
- [Multirate methods](#multirate-methods)
- [Multirate Infinitesimal methods](#infinitesimal-methods)

### Fractional Step integrator

The main function provided is `fractional_step`.

It takes as inputs:

- the operators to use
     	 
	A list of operators (see [Operators](#operators)) in the order you wish to use them.

	Note that the option of a function returning a `Function` will only work with explicit Runge-Kutta methods or analytic subintegration

	If boundary conditions are supplied for explicit, adaptive or EPI methods, the value must be the time derivative of the boundary condition.

	
- delta_t

	the time step to advance by. If using the Irksome solvers, this must be a Constant.

- initial_y

	the initial value(s) to use for y (see [Initial Condition](#initial-condition))
		
- initial_t

	the time value to start at
	
	If using the finite element version, this must be of type Constant
- final_t

	the time value to advance to

- alpha

	the fractional step method to use
	
	this may be a string for one of the pre-defined cases 
	or a user-defined method
	
	For further information, see [Pre-defined fractional step methods](docs/fractional_step.md#defined-fractional-step-methods) or [Defining a new fractional step method](docs/fractional_step.md#defining-a-new-fractional-step-method)

		   
- methods
	
	a dictionary of the methods to use
	
	The key (i, j) will define the method for the jth substep with operator i
	
	The key (i,) will define the default for the ith operator to use if a method is not assigned to a particular substep
	
	The key (0,) will redefine a default to use if no method is assigned to a particular substep or operator
	
	If none of the above are supplied, the default is forward Euler

	There are a variety of built-in options, including [ANALYTIC](docs/fractional_step.md#analytic), [ADAPTIVE](docs/fractional_step.md#adaptive-methods), various [Runge-Kutta methods](docs/fractional_step.md#runge-kutta-methods), various [EPI methods](docs/fractional_step.md#epi-methods), various [EPIRK methods](docs/fractional_step.md#epirk-methods), and methods from [Irksome](docs/fractional_step.md#irksome)

	User defined sub-integrators may also be supplied, of type [Runge-Kutta](docs/fractional_step.md#user-defined-runge-kutta-methods), [Adaptive Runge-Kutta](docs/fractional_step.md#user-defined-adaptive-methods), [EPI](docs/fractional_step.md#user-defined-epi-methods), and [EPIRK](docs/fractional_step.md#user-defined-epirk-methods)

	
Optional Arguments include options for:
- [File Output](docs/fractional_step.md#file-output)
- [Adaptive Fractional Step methods](docs/fractional_step.md#fractional-step-method-control)
- [Dynamic Linearization](docs/fractional_step.md#dynamic-linearization)
- [Boundary Conditions](docs/fractional_step.md#boundary-condition)
- Additional control of the sub-integrators (see the relevant sections of the [Fractional Step documentation](docs/fractional_step.md) or a complete listing of [sub-integrator options](docs/fractional_step.md#options-to-control-sub-integrator-behaviour)).  

For full details, see the complete [fractional step integrator documentation](docs/fractional_step.md)

### Additive Runge Kutta methods

The main function provided is `ark_solve`

It takes as inputs:

- the operators to use
    
	A list of operators (see [Operators](#operators)) in the order you wish to use them.

	They must all either be functions returning arrays or all be Forms

- dt

	the time step to advance by, which may be an integer or a float

- y0

	the initial value(s) to use for y (see [Initial Condition](#initial-condition))

- t0

	the time value to start at
	
	If using the finite element version, this must be of type Constant
- tf

	the time value to advance to

- methods
	a list of `Tableau` or `EmbeddedTableau`, one for each operator.
	
	Note any extra tableaus are ignored.  If there are fewer tableau than operators, no computation is performed.

	See [Additive Runge Kutta method specification](docs/additive_rk.md#method-definition) for further details on how to define a method


Optional Arguments include options for:
- [File Output](docs/additive_rk.md#file-output)
- [Adaptive methods](docs/additive_rk.md#adaptive-methods)
- [Dynamic Linearization](docs/additive_rk.md#dynamic-linearization)
- [Boundary Conditions](docs/additive_rk.md#boundary-condition)
- [Additional control of solver behaviour](docs/additive_rk.md#controlling-solver-behaviour)

For full details, see the [Additive Runge Kutta documentation](docs/additive_rk.md)

### General structure additive runge kutta methods

These are integrated by converting the GARK structure to the ARK structure.

The main function provided is `gark_solve`

It takes as inputs:

- the operators to use
    
	A list of operators (see [Operators](#operators)) in the order you wish to use them.

	They must all either be functions returning arrays or all be Forms

- dt

	the time step to advance by, which may be an integer or a float

- y0

	the initial value(s) to use for y (see [Initial Condition](#initial-condition))

- t0

	the time value to start at
	
	If using the finite element version, this must be of type Constant
- tf

	the time value to advance to

- A
	
	a list of lists containing the arrays a{i, j} defining the method, each as a numpy array

- b
	
	a list containing the vectors b{i} defining the method, each as a numpy array

Optional Arguments include options for:
- [Reordering the method](docs/generalized_additive_rk.md#method-reordering)
- [File Output](docs/generalized_additive_rk.md#file-output)
- [Dynamic Linearization](docs/generalized_additive_rk.md#dynamic-linearization)
- [Boundary Conditions](docs/generalized_additive_rk.md#boundary-condition)
- [Additional control of solver behaviour](docs/generalized_additive_rk.md#controlling-solver-behaviour)

For full details, see the [Generalized Additive Runge Kutta documentation](docs/generalized_additive_rk.md)

### Multirate methods

The multirate methods are integrated by converting to the GARK structure which is then converted to the ARK structure.  The method will be reordered to avoid implicit solves where possible.

The `multirate_solve` function takes as inputs:

- y0

	the initial value(s) to use for y (see [Initial Condition](#initial-condition))

- t0

	the time value to start at
	
	If using the finite element version, this must be of type Constant

- dt

	the time step to advance by, which may be an integer or a float


- tf

	the time value to advance to

- method
	The method to use, which is of type `Multirate`.  See [Multirate Method Definition](docs/multirate.md#method-definition) for further information on defining the method
- M
	the ratio of fast steps to slow steps, which must be an integer
- fs
	the slow operator (see [Operators](#operators))
- ff
	the fast operator (see [Operators](#operators))
	
Optional Arguments include options for:
- [Reordering the method](docs/multirate.md#method-reordering)
- [File Output](docs/multirate.md#file-output)
- [Dynamic Linearization](docs/multirate.md#dynamic-linearization)
- [Boundary Conditions](docs/multirate.md#boundary-condition)
- [Additional control of solver behaviour](docs/multirate.md#controlling-solver-behaviour)

For full details, see the [Multirate documentation](docs/multirate.md)

### Multirate Infinitesimal methods

The `multirate_infinitesimal_solve` function takes as inputs:


- y0

	the initial value(s) to use for y (see [Initial Condition](#initial-condition))

- t0

	the time value to start at
	
	If using the finite element version, this must be of type Constant

- dt

	the time step to advance by, which may be an integer or a float

- tf

	the time value to advance to

- method
	The method to use, which is of type `Multirate_Infinitesimal`.  See [Multirate Infinitesimal Methods](docs/multirate_infinitesimal.md#method-definition) for more details
- fi
	the slow operator (see [Operators](#operators))
	If two slow operators are in use, this is the implicit operator

- ff
	the fast operator (see [Operators](#operators))

- fe (optional)
	the slow explicit operator if two slow operators are in use (see [Operators](#operators))


Optional arguments may be used for:
- [File Output](docs/multirate_infinitesimal.md#file-output)
- [Boundary Conditions](docs/multirate_infinitesimal.md#boundary-condition)
- [Controlling implicit and adaptive stages](docs/multirate_infinitesimal.md#controlling-solver-behaviour)


For full details, see the [Multirate Infinitesimal documentation](docs/multirate_infinitesimal.md)

## Testing Examples

There are a small set of tests available under `examples/tests/`.  There is more complete documentation [here](examples/tests/README.md)

There are also a number of examples demonstrating how to set up problems and use various integrators documented in the [examples directory](examples/Readme.md)
