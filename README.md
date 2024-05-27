# pythOS
	
## Setup

Requirements:

- numpy
- scipy

Optional:

- firedrake - for finite element capabilities
- Irksome - for additional finite element capabilities
- SUNDIALS (v6.7.0) - to use any of the sundials solvers for a subintegration

## Usage

The relevent code must be on the `PYTHONPATH` prior to import.

If using SUNDIALS, prior to use, the wrapper library must be built, and both the wrapper and the SUNDIALS lib folder must be in the `LD_LIBRARY_PATH`

### Fractional Step solver:

The main function provided is `operator_splitting`.

It takes as inputs:

- the operators to use
     	 
	these must be in a list in the order you wish to use them,
	and each must be callable with arguments (t, y) and return a numpy 
	array, unless using the finite element capabilities.
	
	If using firedrake for finite elements, these may be a function 
	callable with arguments (t, y) and returning a `Function` if using 
	only explicit Runge-Kutta methods for the operator. 
	Otherwise, these must be Forms, or tuples containing Forms and boundary
	conditions.  Note that if boundary conditions are supplied for the explicit
	or EPI methods, the value must be the time derivative of the boundary 
	condition.
	
	If using a solver from Irksome, the Form must include the 
	appropriate Dt terms. Otherwise, the Form must not contain any 
	Dt terms.
	
- delta_t

	the time step to advance by, which may be an integer or a float, 
	or a Constant.  If using the Irksome solvers, this must be a Constant.
- initial_y

	the initial value(s) to use for y
	
	this may be a single number or a 1-dimensional numpy array
	
	If using the finite element version, this must be of type Function
	
- initial_t

	the time value to start at
	
	If using the finite element version, this must be of type Constant
- final_t

	the time value to advance to
- alpha

	the operator splitting method to use
	
	this may be a string for one of the pre-defined cases 
	or a user-defined method
	
	The pre-defined 2-splitting methods are
	- 'Godunov': aka Lie-Trotter method, 1st order
	- 'SM2': Strang-Marchuk method, 2nd order, ABBA scheme
	- 'Strang': Strang method, 2nd order, ABAB scheme
	- 'AKOpt22': Auzinger-Ketcheson, optimized 2nd-order, 2-stage method that has minimized local error measure
		https://www.asc.tuwien.ac.at/~winfried/splitting/
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
    - 'AK2s3i-3': Auzinger-Ketcheson 2nd-order 3-stage, ABC scheme v1, https://www.asc.tuwien.ac.at/~winfried/splitting/
    - 'AK2s3ii-3':Auzinger-Ketcheson 2nd-order 3-stage, ABC scheme v2, https://www.asc.tuwien.ac.at/~winfried/splitting/
    - 'AK2s5-3':Auzinger-Ketcheson 2nd-order 5-stage, ABC scheme, https://www.asc.tuwien.ac.at/~winfried/splitting/
	
	Two N-split methods are also defined, which determine 
	the appropriate number of operators based on the input 
	list of operators.
	- 'Godunov-N': Godunov for any number of operators
	- 'Strang-N': Strang for any number of operators

	To define a splitting method, the input is a list of lists.
    The ith entry of the jth list defines the step size for the ith 
	operator in the jth substep.
	
	If a substep need to happen with operators in the reverse order, 
	add an element to the end of the list that is True.
		   
- methods
	
	a dictionary of the methods to use
	
	The key (i, j) will define the method for the jth substep with operator i
	
	The key (i,) will define the default for the ith operator to use if a method is not assigned to a particular substep
	
	The key (0,) will redefine a default to use if no method is assigned to a particular arc or operator
	
	If none of the above are supplied, the default is forward Euler
	
	Built-in options are:
	- ANALYTIC - an analytical solution (or solution using a solver not packaged in pythOS).  To use this, the corresponding operator must be a function that takes as arguments (dt, y) and returns the new solution
	- ADAPTIVE - a solution using an adaptive method to meet some tolerance
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
	- EPI2 - 2nd order EPI method
	- EPI3 - 3rd order EPI method
	- EPI4 - 4th order EPI method
	- EPI5 - 5th order EPI method
	- EPI6 - 6th order EPI method
	- If using Irksome, any of the `ButcherTableau` from Irksome may be used
	
	A user-defined method may also be used by supplying an instance of the `Tableau` class (found in `butcher_tableau.py`) or the `EPIMultistep` class (found in `Epi_multistep.py`) that represent the method.

Optional arguments are:

- b: When using a splitting method which has a parameter, this is the value to use

- fname: a filename to save intermediate results to.  
  
	 When using the finite element capabilities of firedrake this is a 
	 .h5 file accessible through the `CheckpointFile` interface from 
	 `firedrake`. The file has attributes `/times/idx/` that store the time 
	 of each entry and the attribute `/time/last_idx` indicates the last 
	 valid index. The `Function`s are stored using corresponding indices. 
	 Otherwise this is a .csv file containing time in the first entry of each 
	 line, and the solution in the remaining entries

- save\_steps: the number of intermediate steps to save.  
  
	 The default is to save steps at delta\_t time interval if a filename
	 is provided.

- ivp\_methods: a dictionary to control the selection of the ADAPTIVE solver. 
  
	 Key `i` specified the method in use for the ith operator if it is solved 
	 with an ADAPTIVE solver.  The entry has format 
	 (method, relative tolerance, absolute tolerance).
	 Available methods include:
	 - any of the solvers from scipy.integrate.solve\_ivp
	 - Dormand-Prince
	 - Cash-Karp
	 - Fehlberg
	 - Bogacki-Shampine
	 - Heun-Euler
	 - If using SUNDIALS, `CV_ADAMS` and `CV_BDF` for a solver from CVOde,
	 `IDA`, or any of the method strings for MRIStep, ERKStep or ARKStep.
	 
		 If using ARKStep, the method is a tuple specifying (explicit 
		 tableau, explicit tableau), and the provided operator should 
		 also be a tuple of (explicit operator, explicit operator), of 
		 which either (tableau and operator) may be None. 
		 When using MRIStep as an adaptive solver, the provided operator must 
		 have form ((fast implicit, fast explicit), (slow implicit, slow 
		 explicit)) with any of the operators optionally being None. 
		 When using IDA as a solver, the operator must take arguments 
		 (t, y, dydt) and return a residual.
	 - any method not listed, by providing an instance of the `EmbeddedTableau` class found in `butcher_tableau.py`
	 
	When using firedrake, only the adaptive solvers implemented in pythOS 
	are options (or a new method using the `EmbeddedTableau` class)

- epi\_options: a dictionary to control the solver underlying the EPI methods.

	The solver may either be kiops or an ODE integration.  The key i defines 
	the solver to use for the ith operator if it uses an EPI method.  Entries 
	have two forms:
	- ('kiops', tol) to use kiops as the underlying solver
	- (method, (rtol, atol)) - where method is the name of one of the of the integrators from
		scipy.integrate.solve\_ivp, an implemented embedded pair, an instance 
		of the `EmbeddedTableau` class, or `CV_BDF` or `CV_ADAMS` or method 
		from ERKStep or ARKStep (with method specification as with the 
		adaptive options for ARKODE).  This uses an ODE integration as the underlying solver.

- os\_rtol: The relative tolerance for the splitting if using a method with 
  step size control
- os\_atol: The absolute tolerance for the splitting if using a method with
  step size control

- solver_parameters: A dictionary to provide optional arguments to the underlying solvers.

- jacobian: A function to compute a jacobian when using dynamic linearization

	If this option is used, the provided operator functions must take as arguments (t, y, J)

- bc: Any boundary conditions to apply when using firedrake

- stats: A boolean of whether to record stats about accepted or rejected splitting steps when using the splitting methods with adaptive step sizes.

The `Tableau` class is provided in the file `butcher_tableau.py`.  It takes as inputs:

- c : a 1 dimensional numpy array for c
- a : a 2 dimensional numpy array for a
- b : a 1 dimensional numpy array for b

The `EmbeddedTableau` class is provided in the `butcher_tableau.py` file.  It takes as inputs:

- c : a 1 dimensional numpy array for c
- a : a 2 dimensional numpy array for a
- b : a 1 dimensional numpy array for b
- b\_aux : a 1 dimensional numpy array with the embedded method
- order : the lesser of the orders of the two methods, used in step size calculations

### Additive Runge Kutta methods

The `ark_solve` function takes as inputs:

- the operators to use
     	 
	these must be in a list in the order you wish to use them,
	and each must be callable with arguments (t, y) and return a numpy 
	array, unless using the finite element capabilities.
	
	If using firedrake for finite elements, these must be Forms.
- dt

	the time step to advance by, which may be an integer or a float

- y0

	the initial value(s) to use for y
	
	this may be a single number or a 1-dimensional numpy array
	
	If using the finite element version, this must be of type Function

- t0

	the time value to start at
	
	If using the finite element version, this must be of type Constant
- tf

	the time value to advance to

- methods
	a list of `Tableau` or `EmbeddedTableau`, one for each operator.
	
	Note any extra tableaus are ignored.  If there are fewer tableau than operators, no computation is performed.

Optional arguments are:

- fname: a filename to save intermediate results to.  
  
  
	 When using the finite element capabilities of firedrake this is a 
	 .h5 file accessible through the `CheckpointFile` interface from 
	 `firedrake`. The file has attributes `/times/idx/` that store the time 
	 of each entry and the attribute `/time/last_idx` indicates the last 
	 valid index. The `Function`s are stored using corresponding indices. 
	 Otherwise this is a .csv file containing time in the first entry of each 
	 line, and the solution in the remaining entries
 
- save\_steps: the number of intermediate steps to save.  
  
	 The default is to save steps at delta\_t time interval if a filename
	 is provided.

- rtol: The relative tolerance to use if all supplied tableau are `EmbeddedTableau`
- atol: The abolute tolerance to use if all supplied tableau are `EmbeddedTableau`
- bc: Any boundary conditions to apply when using firedrake

- jacobian: A function to compute a jacobian when using dynamic linearization

	If this option is used, the provided operator functions must take as arguments (t, y, J)
	
- solver_parameters: A dictionary to provide optional arguments to the underlying solvers.


### General structure additive runge kutta methods

These are solved by converting the GARK structure to the ARK structure.

The `gark_solve` function takes as inputs:

- the operators to use
     	 
	these must be in a list in the order you wish to use them,
	and each must be callable with arguments (t, y) and return a numpy 
	array, unless using the finite element capabilities.
	
	If using firedrake for finite elements, these must be Forms.
- dt

	the time step to advance by, which may be an integer or a float

- y0

	the initial value(s) to use for y
	
	this may be a single number or a 1-dimensional numpy array
	
	If using the finite element version, this must be of type Function

- initial_t

	the time value to start at
	
	If using the finite element version, this must be of type Constant
- final_t

	the time value to advance to

- A
	
	a list of lists containing the arrays a{i, j} defining the method, each as a numpy array

- b
	
	a list containing the vectors b{i} defining the method, each as a numpy array

Optional arguments are:

- fname: a filename to save intermediate results to.  
  
	 When using the finite element capabilities of firedrake this is a 
	 .h5 file accessible through the `CheckpointFile` interface from 
	 `firedrake`. The file has attributes `/times/idx/` that store the time 
	 of each entry and the attribute `/time/last_idx` indicates the last 
	 valid index. The `Function`s are stored using corresponding indices. 
	 Otherwise this is a .csv file containing time in the first entry of each 
	 line, and the solution in the remaining entries

- save\_steps: the number of intermediate steps to save.  
  
	 The default is to save steps at delta\_t time interval if a filename
	 is provided.

- bc: Any boundary conditions to apply when using firedrake

- jacobian: A function to compute a jacobian when using dynamic linearization

	If this option is used, the provided operator functions must take as arguments (t, y, J)
	
- solver_parameters: A dictionary to provide optional arguments to the underlying solvers.


### Multirate methods

The `multirate_solve` function takes as inputs:

- y0

	the initial value(s) to use for y
	
	this may be a single number or a 1-dimensional numpy array
	
	If using the finite element version, this must be of type Function

- t0

	the time value to start at
	
	If using the finite element version, this must be of type Constant
- dt

	the time step to advance by, which may be an integer or a float

- tf

	the time value to advance to

- method
	
	an instance of the `Multirate` class defining the method to use
	
- fi

	The slow operator.  If the method is IMEX at the slow scale, this is the implicit slow operator.
	
- ff
	
	The fast operator

- fe (optional)

	The explicit slow operator for methods that are IMEX at the slow scale.
	
For each of the operators, they must be callable with arguments (t, y) and return a numpy array, unless using the finite element capabilities.

If using the finite element capabilities, the operators must be Forms.

Optional arguments:

- fname: a filename to save intermediate results to.  
  
	 When using the finite element capabilities of firedrake this is a 
	 .h5 file accessible through the `CheckpointFile` interface from 
	 `firedrake`. The file has attributes `/times/idx/` that store the time 
	 of each entry and the attribute `/time/last_idx` indicates the last 
	 valid index. The `Function`s are stored using corresponding indices. 
	 Otherwise this is a .csv file containing time in the first entry of each 
	 line, and the solution in the remaining entries

- save\_steps: the number of intermediate steps to save.  
  
	 The default is to save steps at delta\_t time interval if a filename
	 is provided.

- ivp\_options: a dictionary of arguments to pass to the adaptive solver.
	
	This includes the tolerances, which have default values of rtol=1e-10 and atol=1e-12

- ivp\_method: the adaptive solver to use.

	The default is RK45 from scipy for the non-finite element version, and 
	Dormand-Prince from the pythOS collection for the finite element version.

- implicit\_solve\_options: A dictionary to provide optional arguments to the underlying solvers.

- bcs: Any Dirichlet boundary conditions to apply in the finite element version.


For convienience, some predefined multirate methods are provided.  These are:

- mri\_kw3 - Knoth and Wolke (1997) order 3
- mri\_erk2a - Sandu (2019) mri-gark-erk22a
- mri\_erk2b - Sandu (2019) mri-gark-erk22b
- mri\_erk3 - Sandu (2019) mri-gark-erk33
- mri\_erk4 - Sandu (2019) mri-gark-erk45a
- mri\_irk2 - Sandu (2019) mri-gark-irk21a
- mri\_esdirk3a - Sandu (2019) mri-gark-esdirk34a
- mri\_sdirk3 - Sandu (2019) mri-gark-sdirk33a
- mri\_imex3 - Chinomona and Reynolds (2021) imex-mri-gark3a
- mri\_imex4 - Chinomona and Reynolds (2021) imex-mri-gark4

The `Multirate` class takes as inputs:

- c - A 1D numpy array
- gamma - A 2D numpy array defining the coupling
- omega (optional) - A 2D numpy array defining the coupling with the second slow operator.
