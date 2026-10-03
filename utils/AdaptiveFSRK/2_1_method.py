from src import fsrk
from sympy import *


# The tableau we expect to get. Taken directly from the Portero paper.
# This is not used at all in calculations, it is only used to compare at the very end.
expected_tableau = fsrk.AdaptiveFSRK2Tableau(
    c=Matrix([1, 1, 4/9, 1/3]),
    A1=Matrix([
        [1, 0, 0, 0],
        [1, 0, 0, 0],
        [-343/180, 0, 47/20, 0],
        [-1592/300, 0, 564/100, 0]
    ]),
    A2=Matrix([
        [0, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 5/9, 0, 0],
        [0, -121/60, 0, 47/20]
    ]),
    baux1=Matrix([1, 0, 0, 0]),
    baux2=Matrix([0, 1, 0, 0]),
    b1=Matrix([1/10, 0, 9/10, 0]),
    b2=Matrix([0, 1/4, 0, 3/4]),
)



# Start with the 2-part fractional implicit Euler scheme
symbolic_tableau = fsrk.AdaptiveFSRK2Tableau(
    c=Matrix([1, 1, Symbol('c_3'), Symbol('c_4')]),
    A1=Matrix([
        [1, 0, 0, 0],
        [1, 0, 0, 0],
        [Symbol('a_3_1'), 0, Symbol('a_3_3'), 0],
        [Symbol('a_4_1'), 0, Symbol('a_4_3'), 0]
    ]),
    A2=Matrix([
        [0, 0, 0, 0],
        [0, 1, 0, 0],
        [0, Symbol('a_3_2'), 0, 0],
        [0, Symbol('a_4_2'), 0, Symbol('a_4_4')]
    ]),
    baux1=Matrix([1, 0, 0, 0]),
    baux2=Matrix([0, 1, 0, 0]),
    b1=Matrix([Symbol('b_1'), 0, Symbol('b_3'), 0]),
    b2=Matrix([0, Symbol('b_2'), 0, Symbol('b_4')]),
)

# Get the order conditions
# TODO: Generate these using something like ketcheson rootedtrees.jl
# Currently these are just hardcoded
conditions = fsrk.get_order_conds(2, symbolic_tableau)

# Impose restrictions from portero paper
# TODO: Generalize! Where do these come from?
e2 = Matrix([0,1,0,0])
e4 = Matrix([0,0,0,1])
restrictions = [
    Eq((symbolic_tableau.A1 * symbolic_tableau.e), symbolic_tableau.c),                 # A1*e = c
    Eq((e2.T * symbolic_tableau.A2 * symbolic_tableau.e)[0,0], symbolic_tableau.c[1]),  # e2^T*A2*e = c1
    Eq((e4.T * symbolic_tableau.A2 * symbolic_tableau.e)[0,0], symbolic_tableau.c[3]),  # e4^T*A2*e = c3 
]

# Combine all equations
conditions = conditions + restrictions

# Get each condition as equal to zero, and filter out any expressions that evaluate to True (these are equations such as 3 = 3 or x = x, not useful for solving the system)
equations = [cond.lhs - cond.rhs for cond in conditions if cond != True]
print("Equations in system to be solved (all =0):")
for eq in equations:
    print(eq)

# Get all symbols in all conditions
unknowns = set()
for eq in equations:
    unknowns.update(eq.free_symbols)

# Remove the symbols we want to remain as free parameters
# TODO: We shouldn't need to ask for specific parameters to be free.
#       The code that comes later needs to not depend on certain free
#       parameters, and when that happens we can remove the block below. 
unknowns.remove(Symbol('b_3'))
unknowns.remove(Symbol('b_4'))
unknowns.remove(Symbol('a_3_3'))
unknowns.remove(Symbol('a_4_3'))
unknowns.remove(Symbol('a_4_4'))

print(f"There are {len(unknowns)} unknowns:")
for un in unknowns:
    print(un)

solutions = solve(equations, unknowns, dict=True)
print("Solutions:")
for symbol in solutions[0]:
    print(f"{symbol} = {solutions[0][symbol]}")

for sol in solutions:
    symbolic_tableau.subs(sol)

# From portero paper: Set a_3_3 = a_4_4 = a to simplify study of A-stability
a_3_3, a_4_4, a = symbols('a_3_3 a_4_4 a')
symbolic_tableau.subs({a_3_3: a, a_4_4: a})

# Get amplification function R & coefficients
R_eq = fsrk.R(symbolic_tableau)
e, deg1, deg2 = fsrk.get_amplification_coeffs(R_eq)

# For this method, a_4_3 is the right variable to zero out the highest degree coefficient
# Will be different for other methods
a_4_3 = Symbol('a_4_3')

# Use the above symbol to obtain nearly l-stable behaviour
l_stability_solution = fsrk.solve_l_stability(e, deg1, deg2, a_4_3)
symbolic_tableau.subs(l_stability_solution[0])

# Recompute R after the l-stability substitution
R_eq_l_stable = fsrk.R(symbolic_tableau)

print("Getting a stability condition, please wait, this may take awhile...")
a_stable = fsrk.get_min_a_stability(R_eq_l_stable, 0.5, grid_size=16)
print(f"a-stability condition: a >= {a_stable}")  # should be close to 2.35 (= 47/20)

c = symbolic_tableau.c
C = diag(*c)
b1, b2 = symbolic_tableau.b1, symbolic_tableau.b2
A1, A2 = symbolic_tableau.A1, symbolic_tableau.A2
e = symbolic_tableau.e

c_sq = Matrix([ci**2 for ci in c])
 
# Terms come from the Portero paper
# TODO: Generalize! Where do these terms come from? Paper does not make this clear.
term1 = 2*((b1.T * c_sq)[0] - Rational(1,3))**2
term2 = 3*((b2.T * c_sq)[0] - Rational(1,3))**2
term3 = ((b1.T * C * A2 * e)[0] - Rational(1,3))**2
term4 = 2*((b1.T * A1 * c)[0] - Rational(1,6))**2
term5 = 3*((b1.T * A2 * c)[0] - Rational(1,6))**2
term6 = 2*((b2.T * A1 * c)[0] - Rational(1,6))**2
term7 = 3*((b2.T * A2 * c)[0] - Rational(1,6))**2
term8 = ((b1.T * A1 * A2 * e)[0] - Rational(1,6))**2
term9 = ((b2.T * A1 * A2 * e)[0] - Rational(1,6))**2

error_coeff = simplify(term1+term2+term3+term4+term5+term6+term7+term8+term9)


b_3 = Symbol('b_3')
a = Symbol('a')

# fix b_4 first
# TODO: Why is b_4 = 3/4 ? Portero paper sets this value, doesn't clearly state how this was found.
substitutions = {Symbol('b_4'): Rational(3, 4)}
error_coeff_final = simplify(error_coeff.subs(substitutions))
print(f"Coefficient of leading term of local error: {error_coeff_final}")

# unconstrained critical point
solutions = solve([diff(error_coeff_final, a), diff(error_coeff_final, b_3)], [a, b_3], dict=True)
a_val = solutions[0][a]

if a_val < a_stable:
    # constraint is binding: fix a = a_stable, re-minimize over b_3 only
    a_val = a_stable
    constrained = error_coeff_final.subs(a, a_stable)
    b_3_val = solve(diff(constrained, b_3), b_3)[0]
else:
    b_3_val = solutions[0][b_3]

print(f"a={a_val}, b_3={b_3_val}")  # expect a=47/20, b_3=9/10

symbolic_tableau.subs(substitutions)
symbolic_tableau.subs({a: a_val, b_3: b_3_val})

print("")
print("EXPECTED TABLEAU:")
expected_tableau.display()
print("")
print("CALCULATED TABLEAU:")
symbolic_tableau.display()
print("")
print("DIFFERENCE (Positive number -> calculated was higher than expected, negative value -> calculated was lower than expected):")
difference_tableau = symbolic_tableau.difference(expected_tableau)
difference_tableau.display()