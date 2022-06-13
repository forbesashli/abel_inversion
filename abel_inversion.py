from pickle import TRUE
import numpy as np
import scipy as sp
import sympy as smp

# test test


class cv:  # cord_variables
    A_VALUE = 3
    X_VALUE = 1
    Y_VALUE = 1


class abel_inversion:
    def abel_inversion():
        x = smp.symbols("x", real=TRUE)
        r = (np.power(cv.A_VALUE, 2) - np.power(cv.Y_VALUE, 2)) ** (1 / 2)
        f = smp.sin(x) ** 3 * smp.exp(-5 * x)
        result = smp.integrate(f, x)

        print(r, result)


abel_inversion.abel_inversion()
