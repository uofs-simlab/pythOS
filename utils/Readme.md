# Utilities for testing various methods

## Work Precision Data collection

The file `work_precision.py` contains a number of functions to collect data for either work-precision diagrams or convergence studies.

The functions within the file are:
- `test_fractional_step`:
    This runs work-precision testing for integrators accessed with the `fractional_step` integrator.  The arguments are:
    - `f_list`: a list of operators
    - `t0`: the initial time
    - `tf`: the final time
    - `y0`: the initial condition
    - `splitting_method`: the desired fractional step method
    - `rk_methods`: the desired integrators for each operator or integration step (as provided to `fractional_step`)
    - `dt0`: the initial delta_t to try for the test
    - `reference_solution`: the reference solution to calculate error
    - `num_time_repeats`: the number of tests to time, of which the minimum will be reported. The default is 0, which will not report any timing
    - `N`: the desired number of data points. The default is 8
    - `err_max`: the maximum error to record. Any tests with larger errors will not be recorded. The default is 0.02
    - `err_min`: the minimum error for the test. Once a solution has a lower error, the test will be stopped and any data collected returned. The default is 2e-12
    - `dt_min`: the minimum delta t to try. Once the testing reaches a delta t smaller, the test will be stopped and any data collected returned.
    - any keyword arguments for `fractional_step`
    The return is a tuple of lists: errors, delta ts, and timings
- `test_additive_rk`:
    This runs work-precision testing for integrators accessed with the `ark_solve` integrator.  The arguments are:
    - `f_list`: a list of operators
    - `t0`: the initial time
    - `tf`: the final time
    - `y0`: the initial condition
    - `tableaus`: the list of tableaus to use
    - `dt0`: the initial delta_t to try for the test
    - `reference_solution`: the reference solution to calculate error
    - `num_time_repeats`: the number of tests to time, of which the minimum will be reported. The default is 0, which will not report any timing
    - `N`: the desired number of data points. The default is 8
    - `err_max`: the maximum error to record. Any tests with larger errors will not be recorded. The default is 0.02
    - `err_min`: the minimum error for the test. Once a solution has a lower error, the test will be stopped and any data collected returned. The default is 2e-12
    - `dt_min`: the minimum delta t to try. Once the testing reaches a delta t smaller, the test will be stopped and any data collected returned.
    - any keyword arguments for `ark_solve`
    The return is a tuple of lists: errors, delta ts, and timings
- `test_gark`:
    This runs work-precision testing for integrators accessed with the `gark_solve` integrator.  The arguments are:
    - `f_list`: a list of operators
    - `t0`: the initial time
    - `tf`: the final time
    - `y0`: the initial condition
    - `A`, `b`: the method specification (as provided to `gark_solve`)
    - `dt0`: the initial delta_t to try for the test
    - `reference_solution`: the reference solution to calculate error
    - `num_time_repeats`: the number of tests to time, of which the minimum will be reported. The default is 0, which will not report any timing
    - `N`: the desired number of data points. The default is 8
    - `err_max`: the maximum error to record. Any tests with larger errors will not be recorded. The default is 0.02
    - `err_min`: the minimum error for the test. Once a solution has a lower error, the test will be stopped and any data collected returned. The default is 2e-12
    - `dt_min`: the minimum delta t to try. Once the testing reaches a delta t smaller, the test will be stopped and any data collected returned.
    - any keyword arguments for `gark_solve`
    The return is a tuple of lists: errors, delta ts, and timings
- `test_multirate`:
    This runs work-precision testing for integrators accessed with the `multirate_solve` integrator.  The arguments are:
    - `fs`: the slow operator
    - `ff`: the fast operator
    - `t0`: the initial time
    - `tf`: the final time
    - `y0`: the initial condition
    - `multirate_method`: the `Multirate` object defining the method
    - `M`: the number of substeps for the fast operator
    - `dt0`: the initial delta_t to try for the test
    - `reference_solution`: the reference solution to calculate error
    - `num_time_repeats`: the number of tests to time, of which the minimum will be reported. The default is 0, which will not report any timing
    - `N`: the desired number of data points. The default is 8
    - `err_max`: the maximum error to record. Any tests with larger errors will not be recorded. The default is 0.02
    - `err_min`: the minimum error for the test. Once a solution has a lower error, the test will be stopped and any data collected returned. The default is 2e-12
    - `dt_min`: the minimum delta t to try. Once the testing reaches a delta t smaller, the test will be stopped and any data collected returned.
    - any keyword arguments for `multirate_solve`
    The return is a tuple of lists: errors, delta ts, and timings
- `test_multirate_infinitesimal`:
    This runs work-precision testing for integrators accessed with the `multirate_solve` integrator.  The arguments are:
    - `fs`: the slow operator
    - `ff`: the fast operator
    - `t0`: the initial time
    - `tf`: the final time
    - `y0`: the initial condition
    - `multirate_method`: the `Multirate_Infinitesimal` object defining the method
    - `dt0`: the initial delta_t to try for the test
    - `reference_solution`: the reference solution to calculate error
    - `fe`: the optional explicit slow operator
    - `num_time_repeats`: the number of tests to time, of which the minimum will be reported. The default is 0, which will not report any timing
    - `N`: the desired number of data points. The default is 8
    - `err_max`: the maximum error to record. Any tests with larger errors will not be recorded. The default is 0.02
    - `err_min`: the minimum error for the test. Once a solution has a lower error, the test will be stopped and any data collected returned. The default is 2e-12
    - `dt_min`: the minimum delta t to try. Once the testing reaches a delta t smaller, the test will be stopped and any data collected returned.
    - any keyword arguments for `multirate_infinitesimal_solve`
    The return is a tuple of lists: errors, delta ts, and timings