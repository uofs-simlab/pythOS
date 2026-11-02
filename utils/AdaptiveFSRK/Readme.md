# Adaptive FSRK Project

Contains tools for generating FSRK methods.

I am currently working towards generalizing the code as much as possible. Currently, if you inspect `2_1_method.py` you will see a lot of code that is arbitrary & specific to this method. The aim to is just give a list of parameters for the method and have it generate according to the specifications.

## Dependencies
- `sympy`, `scipy`, `juliacall` python modules
- `RootedTrees`, `SymPyPythonCall` julia packages
    - These will need to be installed via juliapkg:
    ```shell
    pip install juliapkg

    python3 -m juliapkg add RootedTrees
    python3 -m juliapkg add SymPyPythonCall
    python3 -m juliapkg resolve
    ```

## Checklist
- [x] Be able to replicate the Portero 2(1) method
- [x] Generating order conditions using RootedTrees.jl
- [x] Generating stage order conditions
- [ ] Generalize the a-stability code
- [ ] Generalize the l-stability code
- [ ] Test the code with other Portero FSRK2 methods
- [ ] Potentially generalize to any FSRKm method

# Sources (Not formatted properly yet)
- RootedTrees.jl (Ketcheson)
- Variable step-size fractional step Runge–Kutta methods for time-dependent partial differential equations (Portero)
