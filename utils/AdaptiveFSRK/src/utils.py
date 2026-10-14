def print_equations(eqs, title: str):
    """
    Print a list of sympy equations in a pretty way

    eqs: List of sympy equations
    title: String, a descriptor for this list of equations
    """

    #equations = [eq.lhs - eq.rhs for eq in eqs if eq != True] # Get conditons =0, excluding True statements such as 4=4
    print(title)
    for eq in eqs:
        print(eq)