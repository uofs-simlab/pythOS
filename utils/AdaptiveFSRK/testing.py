from src import fsrk, utils
from sympy import *


# Start with the 2-part fractional implicit Euler scheme
symbolic_tableau = fsrk.AdaptiveFSRKTableau(
    c=Matrix([1, 1, Symbol('c_3', real=True), Symbol('c_4', real=True)]),
    A1=Matrix([
        [1, 0, 0, 0],
        [1, 0, 0, 0],
        [Symbol('a_31', real=True), 0, Symbol('a_33', real=True), 0],
        [Symbol('a_41', real=True), 0, Symbol('a_43', real=True), 0]
    ]),
    A2=Matrix([
        [0, 0, 0, 0],
        [0, 1, 0, 0],
        [0, Symbol('a_32', real=True), 0, 0],
        [0, Symbol('a_42', real=True), 0, Symbol('a_44', real=True)]
    ]),
    baux1=Matrix([1, 0, 0, 0]),
    baux2=Matrix([0, 1, 0, 0]),
    b1=Matrix([Symbol('b_1', real=True), 0, Symbol('b_3', real=True), 0]),
    b2=Matrix([0, Symbol('b_2', real=True), 0, Symbol('b_4', real=True)]),
    order=2
)

utils.print_equations(symbolic_tableau.get_order_conds(), "Order conditions found by RootedTrees.jl:")
