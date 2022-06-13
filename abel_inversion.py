from cgitb import reset
from unittest import result
import numpy as np
import scipy as sp
import sympy as smp
import scipy.integrate as integrate
import scipy.special as special

# test test


class cv:  # cord_variables
    A_VALUE = 2
    X_VALUE = 1.414
    Y_VALUE = 1.414


class A_I:
    def integrand(y, f_r):
        r = smp.symbols("r", real=True)

        return (
            2 * r * f_r / (smp.sqrt(r ** smp.Rational(2, 1) - y ** smp.Rational(2, 1)))
        )

    def abel_inversion(a, y, f_r):
        return smp.integrate(A_I.integrand(y, f_r), (r, y, a)).evalf()

    def test_integrals():

        x = smp.symbols("x", real=True)
        r = (np.power(cv.A_VALUE, 2) - np.power(cv.Y_VALUE, 2)) ** (1 / 2)
        f = smp.sin(x) ** 3 * smp.exp(-5 * x)
        result = integrate.quad(lambda x: special.jv(2.5, x), 0, 4.5)
        print(result)


r = smp.symbols("r", real=True)
f_r = r + 9
a = 2
y = 3
print(A_I.abel_inversion(a, y, f_r))
# A_I.test_integrals()
