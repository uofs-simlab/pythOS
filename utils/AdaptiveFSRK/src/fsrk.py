# fsrk.py

from sympy import *
import numpy as np
from scipy.optimize import minimize
import subprocess
from juliacall import Main as jl

jl.seval("using RootedTrees, SymPyPythonCall")


class AdaptiveFSRKTableau(object):                                           
    def __init__(self, c, A1, A2, baux1, baux2, b1, b2, order):
        """
        All params are sympy Matrices
        c, b1, b2, baux1, baux2 are all 1 dimensional
        A1, A2 are 2 dimensional
        Tableau is expected to be valid on initialization
        """
        self.c, self.A1, self.A2, self.b1, self.b2, self.baux1, self.baux2, self.order = c, A1, A2, b1, b2, baux1, baux2, order

        self.n_stages = shape(c)[0]
        self.e = Matrix([1]*self.n_stages)
    
    def get_order_conds(self):
        """
        Use a Julia script to get order conditions using RootedTrees.jl package
        """
        """
        output = subprocess.run(
            ["julia", "--project=./julia", "./julia/order_conds.jl", srepr(self.A1 + self.A2), srepr(self.baux1 + self.baux2)],
            capture_output=True,
            text=True, # Interpret output as text, not as bytes
            check=True # Throw an exception if julia code fails
        )

        # Parse output
        conditions = []
        for line in output.stdout.splitlines():
            conditions.append(Eq(sympify(line), 0)) # Convert string -> sympy symbols

        return conditions
        """
        jl.seval("""
        using RootedTrees, SymPyPythonCall

        function convert_matrix(A)
            return [SymPyPythonCall.Sym(A[i, j]) for i in 1:size(A, 1), j in 1:size(A, 2)]
        end

        function convert_vector(b)
            return [SymPyPythonCall.Sym(b[i]) for i in 1:length(b)]
        end
        """)
        symbols = {}

        for matrix in (self.A1, self.A2, self.b1, self.b2, self.baux1, self.baux2, self.c):
            for val in matrix:
                if isinstance(val, Symbol):
                    symbols[str(val)] = val

        jl_A1 = jl.convert_matrix(self.A1)
        jl_A2 = jl.convert_matrix(self.A2)
        jl_b1 = jl.convert_vector(self.b1)
        jl_b2 = jl.convert_vector(self.b2)

        ark = jl.AdditiveRungeKuttaMethod([jl_A1, jl_A2], [jl_b1, jl_b2])
        conditions = []

        for p in range(1, self.order+1):
            for t in jl.BicoloredRootedTreeIterator(p):
                residual = jl.residual_order_condition(t, ark)
                eq = Eq(sympify(str(residual), locals=symbols), 0)
                conditions.append(eq)

        return conditions
    
    def subs(self, sub_dict):
        """
        Substitute symbols throughout the whole tableau
        """
        self.c = self.c.subs(sub_dict)
        self.A1 = self.A1.subs(sub_dict)
        self.A2 = self.A2.subs(sub_dict)
        self.b1 = self.b1.subs(sub_dict)
        self.b2 = self.b2.subs(sub_dict)
        self.baux1 = self.baux1.subs(sub_dict)
        self.baux2 = self.baux2.subs(sub_dict)
    
    def difference(self, tableau):
        """
        Compare this tableau with another and return a new tableau containing the differences

        Both tableaux must be the same shape.
        """

        new_c = self.c - tableau.c
        new_A1 = self.A1 - tableau.A1
        new_A2 = self.A2 - tableau.A2
        new_b1 = self.b1 - tableau.b1
        new_b2 = self.b2 - tableau.b2
        new_baux1 = self.baux1 - tableau.baux1
        new_baux2 = self.baux2 - tableau.baux2

        difference_tableau = AdaptiveFSRKTableau(new_c, new_A1, new_A2, new_b1, new_b2, new_baux1, new_baux2, order=self.order)
        return difference_tableau

    
    def display(self):
        row_count = len(self.c.col(0))
        col_count = len(self.A1.row(0))
        length = 0
        for r in range(row_count):
            string = f"{self.c.col(0)[r]:8.4f} |"
            for c in range(col_count):
                string = string + f" {self.A1.row(r)[c]:8.4f}"
            string = string + " |"
            for c in range(col_count):
                string = string + f" {self.A2.row(r)[c]:8.4f}"
            
            print(string)
            length = len(string)

        print("-"*length)

        string = " order 1 |"
        for c in range(col_count):
            string = string + f" {self.baux1.col(0)[c]:8.4f}"
        string = string + " |"
        for c in range(col_count):
            string = string + f" {self.baux2.col(0)[c]:8.4f}"
        print(string)

        string = " order 2 |"
        for c in range(col_count):
            string = string + f" {self.b1.col(0)[c]:8.4f}"
        string = string + " |"
        for c in range(col_count):
            string = string + f" {self.b2.col(0)[c]:8.4f}"
        print(string)


def get_order_conds(p: int, tableau: AdaptiveFSRKTableau):
    """
    Return order conditions equations for adaptive FSRK method with target order p
    Will return a list of sympy equations. If an empty list is returned, it means he function does not support the given order p.
    """

    if p == 2:
        expressions = []

        expr = Eq((tableau.b1.T * tableau.e)[0], 1)
        expressions.append(expr)

        expr = Eq((tableau.b2.T * tableau.e)[0], 1)
        expressions.append(expr)

        expr = Eq((tableau.b1.T * tableau.c)[0], Rational(1, 2)) # Rational(1, 2) == 1 / 2
        expressions.append(expr)

        expr = Eq((tableau.b2.T * tableau.c)[0], Rational(1, 2))
        expressions.append(expr)

        expr = Eq((tableau.b1.T * tableau.A1 * tableau.e)[0], Rational(1, 2))
        expressions.append(expr)

        expr = Eq((tableau.b2.T * tableau.A1 * tableau.e)[0], Rational(1, 2))
        expressions.append(expr)

        expr = Eq((tableau.b1.T * tableau.A2 * tableau.e)[0], Rational(1, 2))
        expressions.append(expr)

        expr = Eq((tableau.b2.T * tableau.A2 * tableau.e)[0], Rational(1, 2))
        expressions.append(expr)
    else:
        return []

    return expressions


def R(tableau: AdaptiveFSRKTableau):
    """
    Get amplification function for an FSRK2 tableau
    Formula below is from Portero paper
    """

    z_1, z_2 = symbols("z_1 z_2")
    I = eye(tableau.n_stages) # Get identity matrix (matrix with 1 along diagonal and 0 elsewhere)

    R_eq = 1 + (z_1 * tableau.b1.T * (I - z_1*tableau.A1 - z_2*tableau.A2).inv() * tableau.e + z_2 * tableau.b2.T * (I - z_1*tableau.A1 - z_2*tableau.A2).inv() * tableau.e)[0]
    
    return cancel(R_eq)


def get_amplification_coeffs(R_eq):
    """
    Given R of some FSRK2 tableau, return coefficients of numerator, e
    Will return dict with keys (i, j), along with deg1, deg2 (So it returns a tuple of length 3)
    """
    z_1, z_2 = symbols("z_1 z_2")
    num, den = fraction(R_eq)

    C = den.subs([(z_1, 0), (z_2, 0)]) # isolates the leftover constant factor
    # Remove constant factor
    num = cancel(num / C)
    den = cancel(den / C)

    num = expand(num)
    num_poly = Poly(num, z_1, z_2)

    deg1 = num_poly.degree(z_1)
    deg2 = num_poly.degree(z_2)

    e = {}
    for i in range(deg1 + 1):
        for j in range(deg2 + 1):
            e[(i, j)] = simplify(num_poly.coeff_monomial(z_1**i * z_2**j))

    return e, deg1, deg2


def solve_l_stability(e, deg1, deg2, target_symbol):
    """
    Zero out the highest-degree numerator term (the (deg1, deg2) corner) to obtain nearly L-stable behaviour

    e: The amplification coefficients
    deg1: highest degree of z_1
    deg2: highest degree of z_2
    target_symbol: Symbol that will be used to zero out the term

    Returns the solution
    """
    high_deg_term = e[(deg1, deg2)]

    sol = solve(Eq(high_deg_term, 0), target_symbol, dict=True)

    return sol


def get_min_a_stability(R_eq, min_val, max_val=10, a_step=0.01, y_range=5, grid_size=4):
    """
    Find the smallest value of 'a' (starting from min_val) where
    |R(i*y1, i*y2, a)| <= 1 for every y1, y2 checked.

    Works by trying increasing values of a, one at a time, and for
    each one checking a grid of (y1, y2) points to see if R ever
    exceeds 1 in magnitude.

    R_eq: the equation for R(z1, z2)
    min_val: smallest value of a tested (will start here and increase)
    max_val: largest value of a tested (will stop here)
    a_step: amount to increase a by on each test
    y_range: the grid of y1,y2 values to test will have each of its two dimensions range from -y_range to +y_range
    grid_size = 4: size of the y1,y2 grid of starting points. If grid_size=4, there will be 16 total starting pairs
                   (minimize() will be used for each pair and check other values beyond these 16)
    """
     
    z_1, z_2 = symbols("z_1 z_2")
    syms = R_eq.free_symbols
    syms.remove(z_1)
    syms.remove(z_2)
    unknown = list(syms)[0]

    # Get R as a regular function (speeds things up)
    R_numeric = lambdify((z_1, z_2, unknown), R_eq, 'numpy')

    # Function that minimize() will use and try to make as small as possible
    def objective(y):
        z1 = 1j * y[0] # recall in python j is the imaginary number
        z2 = 1j * y[1]

        R = R_numeric(z1, z2, value)

        return -abs(R)
    
    # Increase value starting from min_val until we find a stable value
    value = min_val
    while value <= max_val:
        max_magnitude = 0
        starts = np.linspace(-y_range, y_range, grid_size)
        for y1 in starts:
            for y2 in starts:
                y0 = np.array([y1, y2])
                result = minimize(objective, y0)
                max_magnitude = max(max_magnitude, -result.fun) # result.fun is the value that minimize returned

        if max_magnitude <= 1:
            return value
        value += a_step

    return None  # no stable unknown found up to max_val