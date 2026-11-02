def print_equations(eqs, title: str):
    """
    Print a list of sympy equations in a pretty way

    eqs: List of sympy equations
    title: String, a descriptor for this list of equations
    """
    
    print(title)
    for eq in eqs:
        print(eq)