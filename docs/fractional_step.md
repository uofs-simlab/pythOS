# Fractional Step integrator

The fractional step integrator may be used with any number of operators. 

The relevant integrator function is `fractional_step`. It takes some basic arguments to define the problem and method, and optional arguments may be used to control output and sub-integrator behaviour

## Basic Usage

The main function provided is `fractional_step`.

It takes as inputs:

- the operators to use
		 
	A list of operators (see [Operators](../README.md#operators)) in the order you wish to use them.

	Note that the option of a function returning a `Function` will only work with explicit Runge-Kutta methods or analytic subintegration

	If boundary conditions are supplied for explicit, adaptive or EPI methods, the value must be the time derivative of the boundary condition.

	
- delta_t

	the time step to advance by. If using the Irksome solvers, this must be a Constant.

- initial_y

	the initial value(s) to use for y (see [Initial Condition](../README.md#initial-condition))
		
- initial_t

	the time value to start at
	
	If using the finite element version, this must be of type Constant
- final_t

	the time value to advance to

- alpha

	the fractional step method to use
	
	this may be a string for one of the pre-defined cases 
	or a user-defined method
	
	For further information, see [Pre-defined fractional step methods](#defined-fractional-step-methods) or [Defining a new fractional step method](#defining-a-new-fractional-step-method)

		   
- methods
	
	a dictionary of the methods to use
	
	The key (i, j) will define the method for the jth substep with operator i
	
	The key (i,) will define the default for the ith operator to use if a method is not assigned to a particular substep
	
	The key (0,) will redefine a default to use if no method is assigned to a particular substep or operator
	
	If none of the above are supplied, the default is forward Euler

	There are a variety of built-in options, including [ANALYTIC](#analytic), [ADAPTIVE](#adaptive-methods), various [Runge-Kutta methods](#runge-kutta-methods), various [EPI methods](#epi-methods), various [EPIRK methods](#epirk-methods), and methods from [Irksome](#irksome)


## Method Specification
### Defined Fractional step Methods
The pre-defined 2-splitting methods are listed here
- 'Godunov': aka Lie-Trotter method, 1st order
- 'SM2': Strang-Marchuk method, 2nd order, ABBA scheme
- 'Strang': Strang method, 2nd order, ABAB scheme
- 'AKOpt22': Auzinger-Ketcheson, optimized 2nd-order, 2-stage method that has minimized local error measure
    http://www.othmar-koch.org/splitting/
- 'OS22b': can take a parameter b and output a 2nd-order 2-stage method
- 'R3': Ruth, 3rd order
- 'C3': Chambers, 3rd order
- 'S3': Sornborger, 3rd order
- 'AKS3': Auzinger-Ketchenson-Suzuki, 3rd order
- 'AKS3C': Auzinger-Ketchenson-Suzuki with complex coefficients, 3rd order
- 'OS33bv1': can take a parameter b and output a 3nd-order 3-stage method
- 'OS33bv2': can take a parameter b and output a 3nd-order 3-stage method
- 'Y4': Yoshida, 4th order
- 'M4': McLachlan, 4th order
- 'B4': Blanes, 4th order
- 'C4': Chambers, 4th order
- 'AKS3P':

Some 3-spliting methods (ABC schemes):
- 'Godunov-3': Godunov for ABC scheme, 1st order
- 'Strang-3': Strang, 2nd order
- 'Y4-3': Yoshida, 4th order
- 'AK2s3i-3': Auzinger-Ketcheson 2nd-order 3-stage, ABC scheme v1, http://www.othmar-koch.org/splitting/
- 'AK2s3ii-3':Auzinger-Ketcheson 2nd-order 3-stage, ABC scheme v2, http://www.othmar-koch.org/splitting/
- 'AK2s5-3':Auzinger-Ketcheson 2nd-order 5-stage, ABC scheme, http://www.othmar-koch.org/splitting/
	
Two N-split methods are also defined, which determine 
the appropriate number of operators based on the input 
list of operators.
- 'Godunov-N': Godunov for any number of operators
- 'Strang-N': Strang for any number of operators

### Defining a new fractional step method
To define a splitting method, the input is a list of lists.
The ith entry of the jth list defines the step size for the ith 
operator in the jth substep.
	
If a substep need to happen with operators in the reverse order, 
add an element to the end of the list that is True.

## Sub-Integration Method Specification

### Analytic

This option is designed for operators with an analytical solution. 

It may also be used for a sub-integrator that is not packaged in pythOS.

To use this option, the corresponding operator must be a function that takes as arguments (t, dt, y) and returns the new solution

### Adaptive Methods
The key in the `methods` dictionary is simply "ADAPTIVE".  

An optional `ivp_methods` dictionary defines the choice of adaptive solver.  In this case, the selection is made on a per-operator basis.  Key `i` specificies the adaptive method used for the ith operator if an ADAPTIVE solver is used for a sub-integration.  The entry has format (method, relative tolerance, absolute tolerance).
Available methods include:
- any of the solvers from scipy.integrate.solve\_ivp
- Methods built into pythOS:
    - Dormand-Prince
    - Cash-Karp
    - Fehlberg
    - Bogacki-Shampine
    - Heun-Euler
- If the SUNDIALS wrapper and library are available:
    - `CV_ADAMS` or `CV_BDF` from CVOde
    - IDA, which requires an operator function taking arguments (t, y, dydt)
    - ARKStep methods, where the method is specified as (implicit tableau, explicit tableau), and the operator is specified as (implicit operator, explicit operator)
        - for example, ("ARKODE_ARK2_DIRK_3_1_2", "ARKODE_ARK2_ERK_3_1_2").
        - either the implicit tableau or explicit tableau may be None. If one tableau is None, the corresponding operator may also be None.
        - the complete list of tableau is included in [ARKODE documentation](https://sundials.readthedocs.io/en/latest/arkode/Butcher_link.html#)
    - ERKStep methods, which are specified as a single string:
        - for example, "ARKODE_FORWARD_EULER_1_1"
        - the complete list of tableau is included in [ARKODE documentation](https://sundials.readthedocs.io/en/latest/arkode/Butcher_link.html#explicit-butcher-tables)
    - MRIStep methods, where the method is specified as a table ID from [ARKODE](https://sundials.readthedocs.io/en/latest/arkode/Usage/MRIStep/MRIStepCoupling.html#mri-coupling-tables)
        - if using MRIStep, the provided operator must have the form ((fast implicit, fast explicit), (slow implicit, slow explicit))
        - any of the four operators may be None.  In particular, if the slow stage method doesn't have an explicit or implicit component, the corresponding operator should be None
- A user defined method, created as an instance of the `EmbeddedTableau` class

#### Optional method control:
Additional optional arguments may be applied through the `solver_parameters` dictionary, where the options for integrator `i` are supplied through dictionary entry `i` as a dictionary.

The optional arguments available depend on the selection of adaptive integrator:
- Adaptive methods from `scipy.integrate.solve_ivp`
	- optional arguments for `solve_ivp`
- Adaptive methods from SUNDIALS:
	- jac: a jacobian times vector function to use for any method except those from ERKStep and MRIStep
	- max_steps: a maximum number of steps to take in an interval (except for MRIStep)
	- ydot0: an initial value for dydt for IDA
	- id: the vector defining which vector elements are differential (1) or algebraic (0). This is passed to `IDASetId`
	- rtol_slow: the relative tolerance on the slow stepper for MRIStep
	- atol_slow: the absolute tolerator for the slow stepper for MRIStep
	- slow_jac: a jacobian times vector function for the slow operator for MRIStep
	- fast_jac: a jacobian times vector function for the fast operator for MRIStep
	- max_steps_slow: the maximum number of steps to take with the slow operator if using an adaptive slow method for MRIStep
	- max_steps_fast: the maximum number of steps to take with the fast operator per solver interval
	- use_fixed_step: use a fixed step of the size of the solution interval for MRIStep
- Adaptive methods from pythOS built in Embedded Tableaus:
	- bc: a boundary condition to apply if using an operator defined through firedrake
	- any optional arguments (except for a boundary condition) to pass to the firedrake `solve` function if using an operator defined through firedrake
	- any optional arguments to pass to `scipy.optimize.root` (except the `args` argument) if using an implicit integrator and an operator defined by a function returning an array

#### User-Defined Adaptive methods
The `EmbeddedTableau` class is provided in the `butcher_tableau.py` file.  It takes as inputs:

- c : a 1 dimensional numpy array for c
- a : a 2 dimensional numpy array for a
- b : a 1 dimensional numpy array for b
- b\_aux : a 1 dimensional numpy array with the embedded method
- order : the lesser of the orders of the two methods, used in step size calculations


Note that when using Firedrake, only the pythOS adaptive integrators are available, which includes a new user-supplied integrator.

### Runge-Kutta Methods

A number of Runge-Kutta Methods are built in to pythOS
- 'FE' forward Euler
- 'Heun' Heun's method
- 'RK3' the explicit Runge--Kutta of order 3
- 'RK4' the explicit Runge--Kutta of order 4
- 'BE' backward Euler
- 'SD2O2' 2-stage SDIRK method of order 2
- 'SD2O3' 2-stage SDIRK method of order 3
- 'SD3O4' 3-stage SDIRK method of order 4
- 'CN' Crank-Nicolson, also implicit trapezoidal 'IT'
- 'GL2' Gauss-Legendre of order 2, also implicit midpoint 'IM'
- 'GL4' Gauss-Legendre of order 4
- 'GL6' Gauss-Legendre of order 6
- 'SSP(5,4)'
- 'SSPRK3' Strong stability preserving RK3
- 'SDAstable' 2-stage, second order SDIRK method (Pareschi and Russo with x=1/2)
- 'SDLstable' 2 stage, second order SDIRK method (Pareschi and Russo with x=(2 + sqrt(2))/2)
- 'SD3O3Lstable' 3 stage SDIRK method of order 3
- 'SD5O4' 5-stage SDIRK method of order 4

In addition, a user-supplied method may be supplied by supplying an instance of the `Tableau` class

#### Optional method control:
Additional optional arguments may be applied through the `solver_parameters` dictionary, where the options for integrator `i` are supplied through dictionary entry `i` as a dictionary.

- bc: a boundary condition to apply if using an operator defined through firedrake
- any optional arguments (except for a boundary condition) to pass to the firedrake `solve` function if using an operator defined through firedrake
- any optional arguments to pass to `scipy.optimize.root` (except the `args` argument) if using an implicit integrator and an operator defined by a function returning an array

#### User defined Runge-Kutta methods

The `Tableau` class is provided in the file `butcher_tableau.py`.  It takes as inputs:

- c : a 1 dimensional numpy array for c
- a : a 2 dimensional numpy array for a
- b : a 1 dimensional numpy array for b

Depending on the structure of `a`, pythOS will use an explicit, diagonally implicit, or fully implicit Runge Kutta integrator

### EPI methods

There is the option of using exponential propagation iterative methods (EPI).  Because they are multistep methods, by default, they substep on each use, and reset the history after each solve interval.  
If you wish to avoid resetting the history, the optional `monolithic` argument can be set to true (through the `solver_parameters` dictionary).  This is only recommended if the use is truly monolithic (i.e. no splitting is applied, and the step size is constant). To build the history, the initial steps are solved with a highly accurate adaptive method (relative tolerance=1e-10, absolute tolerance=1e-12).
The number of substeps can be controlled with the `n_steps` optional argument (through the `solver_parameters` dictionary). The default number is 20.

Built in EPI methods are:
- EPI2 - 2nd order EPI (multistep) method 
- EPI3 - 3rd order EPI (multistep) method
- EPI4 - 4th order EPI (multistep) method
- EPI5 - 5th order EPI (multistep) method
- EPI6 - 6th order EPI (multistep) method


Additionally, the underlying solver may be selected by using the optional `epi_options` dictionary. The solver may either be kiops or an ODE integration.  The key i defines 
the solver to use for the ith operator if it uses an EPI method.  Entries 
have two forms:
- ('kiops', tol) to use kiops as the underlying solver
- (method, (rtol, atol)) - where method is the name of one of the of the integrators listed in [Adaptive Methods](#adaptive-methods), excluding the MRI methods.  If an ARKStep tuple of methods is specified, the implicit solver is used.

#### Optional method control:
Additional optional arguments may be applied through the `solver_parameters` dictionary, where the options for integrator `i` are supplied through dictionary entry `i` as a dictionary.
- n_steps: the number of sub-steps to take in each step
- monolithic: if True, the solver will not reset the history unless the size of the time step changes
- bc: an optional boundary condition to apply if the operator definition is defined with Firedrake
- kiops options:
	- m_init: (kiops only) - initial m (Krylov size)
	- mmin: (kiops only) - minimum m (Krylov size)
	- mmax: (kiops only) - maximum m (Krylov size)
	- iop: (kiops only) - length of the incomplete orthogonalization procedure
	- optional arguments for firedrake `solve` function if using kiops and a firedrake operator
- exode options:
	- any optional arguments available for the underlying adaptive integrator detailed above

#### User Defined EPI methods
A new method may be created by creating an instance of the `EpiMultistep` class.  This takes as an argument an array `A`, representing the coefficients of the method. Note that this doesn't include the zeroed first row.  A single instance of the `EpiMultistep` class should not be used for two operators, because it also stores the previously computed values necessary for advancing the solution.


### EPIRK methods

There is the option of using exponential propagation iterative methods of type Runge-Kutta (EPIRK).  
A number of EPIRK methods are built into pythOS:
- EPIRK2 - 2nd order EPIRK method
- EPIRK3 - 3rd order EPIRK method
- EPIRK4 - 4th order EPIRK method
- EPIRK4s3 - 4th order EPIRK method
- EPIRK5s3 - 5th order EPIRK method (with 3 stages)
- EPIRK5s4 - 5th order EPIRK method (with 4 stages)

Additionally, the underlying solver may be selected by using the optional `epi_options` dictionary. The solver may either be kiops or an ODE integration.  The key i defines 
the solver to use for the ith operator if it uses an EPI method.  Entries 
have two forms:
- ('kiops', tol) to use kiops as the underlying solver
- (method, (rtol, atol)) - where method is the name of one of the of the integrators listed in [Adaptive Methods](#adaptive-methods), excluding the MRI methods.  If an ARKStep tuple of methods is specified, the implicit solver is used.

#### Optional method control:
Additional optional arguments may be applied through the `solver_parameters` dictionary, where the options for integrator `i` are supplied through dictionary entry `i` as a dictionary.
- bc: an optional boundary condition to apply if the operator definition is defined with Firedrake
- kiops options:
	- m_init: (kiops only) - initial m (Krylov size)
	- mmin: (kiops only) - minimum m (Krylov size)
	- mmax: (kiops only) - maximum m (Krylov size)
	- iop: (kiops only) - length of the incomplete orthogonalization procedure
	- optional arguments for firedrake `solve` function if using kiops and a firedrake operator
- exode options:
	- any optional arguments available for the underlying adaptive integrator detailed above


#### User defined EPIRK methods
A new method may be created by creating an instance of the `EpiRKMethod` class.  This takes as an argument arrays `p`, `g`, and `ab`.  Note that `ab` is an array containing the coefficients of a, and b in the last row.

### Irksome

Methods from Irksome may be used when the problem being solved is expressed with finite element discretization.
If you wish to use Irksome, the operator definition must include all the relevant Dt terms.
In this case, the method specification is the `ButcherTableau` object from Irksome.
In addition, when Irksome is in use, the delta t, and initial t arguments must be `Constant`s (or `MeshConstant.Constant`s)

If a split method is desired, or any other optional arguments to `TimeStepper` are desired, these can be supplied on a per-operator basis to the `solver_parameters` dictionary.  The entry should be a dictionary containing the optional arguments. For further information on the optional arguments, see the [Irksome documentation](https://www.firedrakeproject.org/Irksome/index.html)

## Optional Arguments

Some optional arguments are specific to the selection of sub-integrator, which are detailed in the sections above

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
	 Note that if an adaptive fractional step method is used, the saved steps may not be evenly spaced

### Fractional Step method control
- b: When using a splitting method which has a parameter, this is the value to use
- os\_rtol: The relative tolerance for the splitting if using a method with 
  step size control
- os\_atol: The absolute tolerance for the splitting if using a method with
  step size control
- stats: A boolean of whether to record stats about accepted or rejected splitting steps when using the splitting methods with adaptive step sizes.

### Dynamic Linearization

- jacobian: A function to compute a jacobian when using dynamic linearization

	If this option is used, the provided operator functions must take as arguments (t, y, J)

When this is used, the jacobian function is called at the start of each step, and held constant for each substep within.  This jacobian is passed to each operator function, which may be uses to compute the values of each operator

### Boundary Condition
- bc: Any boundary conditions to apply when using firedrake

This boundary condition will be enforced after each time step if it is supplied.

### Options to control sub-integrator behaviour
- ivp\_methods: a dictionary to control the selection of the ADAPTIVE solver. 
  
	 Key `i` specified the method in use for the ith operator if it is solved 
	 with an ADAPTIVE solver.  See [Adaptive methods](#adaptive-methods) for more details about the contents
- epi\_options: a dictionary to control the solver underlying the EPI or EPIRK methods.
	 Key `i` specifies the solver selection, which is either kiops or an ODE integration.  See [EPI methods](#epi-methods) or [EPIRK methods](#epirk-methods) for further details

- solver_parameters: A dictionary to provide optional arguments to the underlying solvers.
	 Key `i` specifies the optional arguments for the solvers underlying the sub-integrator for operator `i`.
	 
	 The available optional arguments depend on the sub-integrator and operator definition.

	 - Runge Kutta methods:
	 	- bc: a boundary condition to apply if using an operator defined through firedrake
		- any optional arguments (except for a boundary condition) to pass to the firedrake `solve` function if using an operator defined through firedrake
		- any optional arguments to pass to `scipy.optimize.root` (except the `args` argument) if using an implicit integrator and an operator defined by a function returning an array
	 - Adaptive methods from pythOS built in Embedded Tableaus:
	 	The same optional arguments as the Runge Kutta methods. They are currently all explicit, but that is not necessarily required
	 - Adaptive methods from `scipy.integrate.solve_ivp`
	 	- optional arguments for `solve_ivp`
	 - Adaptive methods from SUNDIALS:
	 	- jac: a jacobian times vector function to use for any method except those from ERKStep and MRIStep
		- max_steps: a maximum number of steps to take in an interval (except for MRIStep)
		- ydot0: an initial value for dydt for IDA
		- id: the vector defining which vector elements are differential (1) or algebraic (0). This is passed to `IDASetId`
		- rtol_slow: the relative tolerance on the slow stepper for MRIStep
		- atol_slow: the absolute tolerator for the slow stepper for MRIStep
		- slow_jac: a jacobian times vector function for the slow operator for MRIStep
		- fast_jac: a jacobian times vector function for the fast operator for MRIStep
		- max_steps_slow: the maximum number of steps to take with the slow operator if using an adaptive slow method for MRIStep
		- max_steps_fast: the maximum number of steps to take with the fast operator per solver interval
		- use_fixed_step: use a fixed step of the size of the solution interval for MRIStep
	 - EPI methods: 
	 	- n_steps: the number of sub-steps to take in each step
		- monolithic: if True, the solver will not reset the history unless the size of the time step changes
		- bc: an optional boundary condition to apply if the operator definition is defined with Firedrake
		- m_init: (kiops only) - initial m (Krylov size)
		- mmin: (kiops only) - minimum m (Krylov size)
		- mmax: (kiops only) - maximum m (Krylov size)
		- iop: (kiops only) - length of the incomplete orthogonalization procedure
		- optional arguments for firedrake `solve` function if using kiops and a firedrake operator
		- optional arguments for `scipy.integrate.solve_ivp` if using exode solver and a scipy integrator
		- optional arguments detailed above for adaptive solvers from SUNDIALS (for exode with a Sundials solver) or from pythOS (for adaptive solvers built in to pythOS)
	 - EPIRK methods:
		- bc: an optional boundary condition to apply if the operator definition is defined with Firedrake
		- m_init: (kiops only) - initial m (Krylov size)
		- mmin: (kiops only) - minimum m (Krylov size)
		- mmax: (kiops only) - maximum m (Krylov size)
		- iop: (kiops only) - length of the incomplete orthogonalization procedure
		- optional arguments for firedrake `solve` function if using kiops and a firedrake operator
		- optional arguments for `scipy.integrate.solve_ivp` if using exode solver and a scipy integrator
		- optional arguments detailed above for adaptive solvers from SUNDIALS (for exode with a Sundials solver) or from pythOS (for adaptive solvers built in to pythOS)
	 - Irksome methods:
		- optional arguments to pass to `TimeStepper` from Irksome if the operator and sub-integrator are from Irksome



