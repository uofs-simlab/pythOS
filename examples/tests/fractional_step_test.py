from brusselator_ode_problem import *
import fractional_step as fs
import sys
from testing_utils import plot, output

verbose = "-v" in sys.argv or "--verbose" in sys.argv
enable_sundials = "-s" in sys.argv or "--sundials" in sys.argv

if len(sys.argv) == 2 and (sys.argv[1] == "-h" or sys.argv[1] == '--help'):
    print("Usage:")
    print("python3 multirate_test.py [-p|--plot] [-v|--verbose] [-s|--sundials] [-h|--help]")
    print("""Options:
    -p | --plot\tProduce plot of the results
    -v | --verbose\tPrint the solution at tf for all methods tested
    -s | --sundials\tTry the sundials adaptive solver
    -h | --help\tDisplay this help message and exit (if it is the only command line flag)""")
    exit()


solution = fs.fractional_step([lambda t, y: f1(t,y) + f2(t,y)], .1, y0, 0, tf, 'Godunov', {(0,): 'ADAPTIVE'}, fname='adaptive.csv')
if verbose:
    print("{:<40} {}".format("adaptive solution", solution))

## Basic solvers
result = fs.fractional_step([f1,f2], 0.1, y0, 0, tf, 'R3', methods={(2,): "SD2O3",(1,):"RK3"}, fname='split.csv')

output(verbose, 0.1, result, solution, "basic split")

result = fs.fractional_step([f1, f2], 0.05, y0, 0, tf, 'Strang', methods={(2,): "CN", (1,): "Heun"},fname='split_fully_implicit.csv')
output(verbose, 0.1, result, solution, "using fully implicit solver")

## Adaptive solvers
result = fs.fractional_step([f1,f2], .1, y0, 0, tf, 'Strang', {(1,): 'ADAPTIVE', (2,): "SD2O2"}, fname='adaptive_split.csv')

output(verbose, 0.05, result, solution, "using adaptive solver")

result = fs.fractional_step([f1,f2], .1, y0, 0, tf, 'Strang', {(1,): 'ADAPTIVE', (2,): "SD2O2"}, ivp_methods={1: ("Cash-Karp", 1e-6, 1e-8)}, fname='adaptive_Cash-Karp_split.csv')
output(verbose, 0.05, result, solution, "using Cash-Karp solver")

if enable_sundials:
    result = fs.fractional_step([f1, f2], .1, y0, 0, tf, "Strang", {(1,): "ADAPTIVE", (2,): "SD2O2"}, ivp_methods={1: ("CV_BDF", 1e-6, 1e-8)}, fname='adaptive_CVODE_split.csv')
    output(verbose, 0.05, result, solution, "using CVODE solver")

## EPI solvers
result = fs.fractional_step([f1,f2], 0.1, y0, 0, tf, 'Strang', {(2,): "EPI2",(1,):"RK4"}, epi_options={2: ('kiops', 1e-10)}, fname='epi.csv')
output(verbose, 0.1, result, solution, "using EPI (kiops) solver")

result = fs.fractional_step([f1,f2], 0.1, y0, 0, tf, 'Strang', {(2,): "EPI2",(1,): "RK4"},epi_options={2: ('RK45', (1e-6,1e-8))},fname='epi_RK45.csv')
output(verbose, 0.05, result, solution, "using EPI (RK45) solver")

result = fs.fractional_step([f1,f2], 0.1, y0, 0, tf, 'Strang', {(2,): "EPI2",(1,): "RK4"},epi_options={2: ('Cash-Karp', (1e-6,1e-8))},fname='epi_Cash-Karp.csv')
output(verbose, 0.05, result, solution, "using EPI (Cash-Karp) solver")

if enable_sundials:
    result = fs.fractional_step([f1, f2], 0.1, y0, 0, tf, 'Strang', {(2,): "EPI2", (1,): "RK4"}, epi_options={2: ('ARKODE_FEHLBERG_6_4_5', (1e-6,1e-8))},fname='epi_ARKODE.csv')
    output(verbose, 0.05, result, solution, "using EPI (ARKODE) solver")

## Basic (complex) solvers
result = fs.fractional_step([f1,f2], 0.1, y0, 0, tf, 'C4', methods={(2,): "SD2O3",(1,):"RK3"}, fname='split_complex.csv')

output(verbose, 0.1, result, solution, "basic (complex) split")

result = fs.fractional_step([f1, f2], 0.05, y0, 0, tf, 'C4', methods={(2,): "CN", (1,): "Heun"},fname='split_fully_implicit_complex.csv')
output(verbose, 0.1, result, solution, "using fully implicit solver")

## Adaptive solvers
result = fs.fractional_step([f1,f2], .1, y0, 0, tf, 'C4', {(1,): 'ADAPTIVE', (2,): "SD2O2"}, ivp_methods={1: ("Cash-Karp", 1e-6, 1e-8)}, fname='adaptive_Cash-Karp_split_complex.csv')
output(verbose, 0.05, result, solution, "using Cash-Karp solver")

## EPI solvers
result = fs.fractional_step([f1,f2], 0.1, y0, 0, tf, 'C4', {(2,): "EPI2",(1,):"RK4"}, epi_options={2: ('kiops', 1e-10)}, fname='epi_complex.csv')
output(verbose, 0.1, result, solution, "using EPI (kiops) solver")

result = fs.fractional_step([f1,f2], 0.1, y0, 0, tf, 'C4', {(2,): "EPI2",(1,): "RK4"},epi_options={2: ('RK45', (1e-6,1e-8))},fname='epi_RK45_complex.csv')
output(verbose, 0.05, result, solution, "using EPI (RK45) solver")

result = fs.fractional_step([f1,f2], 0.1, y0, 0, tf, 'C4', {(2,): "EPI2",(1,): "RK4"},epi_options={2: ('Cash-Karp', (1e-6,1e-8))},fname='epi_Cash-Karp_complex.csv')
output(verbose, 0.05, result, solution, "using EPI (Cash-Karp) solver")

    
if '-p' in sys.argv or '--plot' in sys.argv:
    basic_labels = {
        'split.csv': "basic split",
        'split_fully_implicit.csv': "split with fully implicit solver",
        'split_complex.csv': "basic complex split",
        'split_fully_implicit_complex.csv': "complex split with fully implicit solver"
    }
    adaptive_labels = {
        'adaptive_split.csv': "split with adaptive solver",
        'adaptive_Cash-Karp_split.csv': "split with Cash-Karp solver",
        'adaptive_Cash-Karp_split_complex.csv': "complex split with Cash-Karp solver"
    }
    epi_labels = {
        'epi.csv': "split with EPI (kiops) solver",
        'epi_RK45.csv': "split with EPI (RK45) solver",
        'epi_Cash-Karp.csv': "split with EPI (Cash-Karp) solver",
        'epi_complex.csv': "complex split with EPI (kiops) solver",
        'epi_RK45_complex.csv': "complex split with EPI (RK45) solver",
        'epi_Cash-Karp_complex.csv': "complex split with EPI (Cash-Karp) solver"
    }
    
    if enable_sundials:
        adaptive_labels["adaptive_CVODE_split.csv"] = 'split with CVODE solver'
        epi_labels['epi_ARKODE.csv'] = 'split with EPI (ARKODE) solver'

    plot('adaptive.csv', basic_labels, title='Basic Solvers',show=False)
    plot('adaptive.csv', adaptive_labels, title='Splits using Adaptive Solvers',show=False)
    plot('adaptive.csv', epi_labels, title='Splits using EPI solvers')
