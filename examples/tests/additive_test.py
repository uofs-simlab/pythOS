import additive_rk as ark
import gark_methods as gark
from butcher_tableau import Tableau, EmbeddedTableau
import numpy as np
from brusselator_ode_problem import *
from testing_utils import plot, output
import sys
import fractional_step as fs

verbose = "-v" in sys.argv or "--verbose" in sys.argv
if len(sys.argv) == 2 and (sys.argv[1] == "-h" or sys.argv[1] == '--help'):
    print("Usage:")
    print("python3 multirate_test.py [-p|--plot] [-v|--verbose] [-h|--help]")
    print("""Options:
    -p | --plot\tProduce plot of the results
    -v | --verbose\tPrint the solution at tf for all methods tested
    -h | --help\tDisplay this help message and exit (if it is the only command line flag)""")
    exit()

tableau_3_e = Tableau(
    np.array([0, 1767732205903/2027836641118, 3/5, 1]),
    np.array([[0,0,0,0],
              [1767732205903/2027836641118,0,0,0],
              [5535828885825/10492691773637,788022342437/10882634858940,0,0],
              [6485989280629/16251701735622, -4246266847089/9704473918619,10755448449292/10357097424841,0]]),
    np.array([1471266399579/7840856788654, -4482444167858/7529755066697, 11266239266428/11593286722821, 1767732205903/4055673282236]))
tableau_3_i = Tableau(
    np.array([0, 1767732205903/2027836641118, 3/5, 1]),
    np.array([[0,0,0,0],
              [1767732205903/4055673282236,1767732205903/4055673282236,0,0],
              [2746238789719/10658868560708, -640167445237/6845629431997, 1767732205903/4055673282236, 0],
              [1471266399579/7840856788654, -4482444167858/7529755066697, 11266239266428/11593286722821, 1767732205903/4055673282236]]),
    np.array([1471266399579/7840856788654, -4482444167858/7529755066697, 11266239266428/11593286722821, 1767732205903/4055673282236]))

solution = fs.fractional_step([lambda t, y: f1(t,y) + f2(t,y)], .1, y0, 0, tf, 'Godunov', {(0,): 'ADAPTIVE'}, fname='adaptive.csv')
if verbose:
    print("{:<40} {}".format("adaptive solution", solution))

result = ark.ark_solve([f1,f2], 0.1, y0, 0, tf, [tableau_3_e,tableau_3_i],fname='ark_solve.csv')
output(verbose, 0.02, result, solution, "Additive RK")


tableau_e = EmbeddedTableau(
    np.array([0, 1/2, 83/250, 31/50, 17/20, 1]),
    np.array([[0,0,0,0,0,0],
              [1/2,0,0,0,0,0],
              [13861/62500, 6889/62500,0,0,0,0],
              [-116923316275/2393684061468, -2731218467317/15368042101831, 9408046702089/11113171139209, 0, 0, 0],
              [-451086348788/2902428689909, -2682348792572/7519795681897, 12662868775082/11960479115383, 3355817975965/11060851509271, 0, 0],
              [647845179188/3216320057751, 73281519250/8382639484533, 552539513391/3454668386233, 3354512671639/8306763924573, 4040/17871 , 0]]),
    np.array([82889/524892, 0, 15625/83664, 69875/102672, -2260/8211, 1/4]),
    np.array([4586570599.0/29645900160.0, 0, 178811875.0/945068544.0, 814220225.0/1159782912.0, -3700637.0/11593932.0, 61727.0/225920.0]), 3)
tableau_i = EmbeddedTableau(
    np.array([0, 1/2, 83/250, 31/50, 17/20, 1]),
    np.array([[0,0,0,0,0,0],
              [1/4, 1/4, 0, 0 ,0 ,0],
              [8611/62500, -1743/31250 ,1/4 ,0 ,0 ,0],
              [5012029/34652500, -654441/2922500, 174375/388108, 1/4, 0, 0],
              [15267082809/155376265600 ,-71443401/120774400 ,730878875/902184768 ,2285395/8070912 ,1/4 ,0],
              [82889/524892, 0, 15625/83664, 69875/102672, -2260/8211, 1/4]]),
    np.array([82889/524892, 0, 15625/83664, 69875/102672, -2260/8211, 1/4]),
    np.array([4586570599.0/29645900160.0, 0, 178811875.0/945068544.0, 814220225.0/1159782912.0, -3700637.0/11593932.0, 61727.0/225920.0]), 3)

result = ark.ark_solve([f1,f2], 0.2, y0, 0, tf, [tableau_e,tableau_i],fname='ark_solve_adaptive.csv',rtol=1e-6,atol=1e-8)
output(verbose, 0.02, result, solution, "Adaptive additive RK")

beta = 0.5
Aee = np.array([[0, 0, 0], [1/2, 0, 0], [1-beta, beta, 0]])
Aei = np.array([[0, 0], [1/2, 0], [1/2, 1/2]])
Aie = np.array([[1/4, 0, 0], [1/4, 1/2, 0]])
Aii = np.array([[1/4, 0], [1/2, 1/4]])

bE = np.array([1/4, 1/2, 1/4])
bI = np.array([1/2, 1/2])

A = [[Aee, Aei], [Aie, Aii]]
b = [bE, bI]
result = gark.gark_solve([f1,f2],0.1,y0,0,tf,A,b, fname='gark_solve.csv')
output(verbose,0.1,result,solution, "Generalized additive RK")

if '-p' in sys.argv or '--plot' in sys.argv:
    labels={'ark_solve.csv': "Additive RK",
            'ark_solve_adaptive.csv': "Adaptive additive RK",
            'gark_solve.csv': "Generalized additive RK"}
    plot('adaptive.csv', labels)
