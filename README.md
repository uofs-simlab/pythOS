# pythOS

pythOS is a library for solving differential equations through operator splitting.

## Setup

Requirements:

- numpy
- scipy

Optional:

- firedrake - for finite element capabilities

## Usage

The relevent code must be on the `PYTHONPATH` prior to import.

### Fractional Step solver:

The main function provided is `operator_splitting`.

It takes as inputs:

- the operators to use
     	 
	these must be in a list in the order you wish to use them,
	and each must be callable with arguments (t, y) and return a numpy 
	array, unless using the finite element capabilities.
	
	If using firedrake for finite elements, these must be Forms. 
		
- delta_t

	the time step to advance by, which may be an integer or a float.
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
		https://www.asc.tuwien.ac.at/auzinger/splitting/index.php
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
	- 'AKS3P': Auzinger 3-stage 3rd order complex scheme

    Some 3-spliting methods (ABC schemes):
    - 'Godunov-3': Godunov for ABC scheme, 1st order
    - 'Strang-3': Strang, 2nd order
    - 'Y4-3': Yoshida, 4th order
    - 'AK2s3i-3': Auzinger-Ketcheson 2nd-order 3-stage, ABC scheme v1, https://www.asc.tuwien.ac.at/auzinger/splitting/index.php
    - 'AK2s3ii-3':Auzinger-Ketcheson 2nd-order 3-stage, ABC scheme v2, https://www.asc.tuwien.ac.at/auzinger/splitting/index.php
    - 'AK2s5-3':Auzinger-Ketcheson 2nd-order 5-stage, ABC scheme, https://www.asc.tuwien.ac.at/auzinger/splitting/index.php
	
	Two N-split methods are also defined, which determine 
	the appropriate number of operators based on the input 
	list of operators.
	- 'Godunov-N': Godunov for any number of operators
	- 'Strang-N': Strang for any number of operators

	To define a splitting method, the input is a list of lists.
    The ith entry of the jth list defines the step size for the ith 
	operator in the jth substep.
			   
- methods
	
	a dictionary of the methods to use
	
	The key (i, j) will define the method for the jth substep with operator i
	
	The key (i,) will define the default for the ith operator to use if a method is not assigned to a particular substep
	
	The key (0,) will redefine a default to use if no method is assigned to a particular arc or operator
	
	If none of the above are supplied, the default is forward Euler
	
	Built-in options are:
	- ANALYTIC - an analytical solution (or solution using a solver not packaged in pythOS).  To use this, the corresponding operator must be a function that takes as arguments (dt, y) and returns the new solution
	- EXACT - a numerically exact solution.
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
	
	A user-defined method may also be used by supplying an instance of the `Tableau` class (found in `butcher_tableau.py`).

Optional arguments are:

- b: When using a splitting method which has a parameter, this is the value to use

- fname: a filename to save intermediate results to.  
  
	  When using the finite element capabilities of firedrake this is a 
	  .h5 file accessible through the `CheckpointFile` interface from 
	  firedrake. Otherwise this is a .csv file containing time in the 
	  first entry of each line, and the solution in the remaining entries

- save\_steps: the number of intermediate steps to save.  
  
	  The default is to save steps at delta\_t time interval if a filename
	  is provided.

- ivp\_methods: a dictionary to control the selection of the EXACT solver. 
  
	 Key `i` specified the method in use for the ith operator if it is solved 
	 with an EXACT solver.  The entry has format 
	 (method, relative tolerance, absolute tolerance).
	 Available methods include:
	 - any of the solvers from scipy.integrate.solve\_ivp
	 - Dormand-Prince
	 - Cash-Karp
	 - Fehlberg
	 - Bogacki-Shampine
	 - Heun-Euler
	 - any method not listed, by providing an instance of the `EmbeddedTableau` class found in `butcher_tableau.py`
	 
	  When using firedrake, only the exact solvers implemented in pythOS 
	  are options (or a new method using the `EmbeddedTableau` class)

- solver_parameters: A dictionary to provide optional arguments to the underlying solvers.

- bc: Any boundary conditions to apply when using firedrake

