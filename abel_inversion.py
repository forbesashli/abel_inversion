import numpy as np
import scipy as sp
import sympy as smp
from sympy import I
import scipy.integrate as integrate
import scipy.special as special
import yaml

# test test


def readInConfigFile(ymlFile):
    with open(ymlFile, "r") as stream:
        try:
            configDict = yaml.safe_load(stream)
            return configDict
        except yaml.YAMLError as exc:
            print(exc)
            return None


CONFIG_DICT = readInConfigFile("config.yaml")


class cv:  # cord_variables
    A_VALUE = CONFIG_DICT["INTEGRAL_VARIABLES"]["A_VAL"]
    R_VALUE = CONFIG_DICT["INTEGRAL_VARIABLES"]["R_VAL"]
    X_VALUE = CONFIG_DICT["INTEGRAL_VARIABLES"]["X_VAL"]
    Y_VALUE = CONFIG_DICT["INTEGRAL_VARIABLES"]["Y_VAL"]


class A_T:  # abel transform
    def integrand(y, f_r):
        r = smp.symbols("r", real=True)

        return (
            2 * r * f_r / (smp.sqrt(r ** smp.Rational(2, 1) - y ** smp.Rational(2, 1)))
        )

    def abel_transform(a, y, f_r):
        r = smp.symbols("r", real=True)

        return smp.integrate(A_T.integrand(y, f_r), (r, y, a)).evalf()

    def test_integrals():

        x = smp.symbols("x", real=True)
        r = (np.power(cv.A_VALUE, 2) - np.power(cv.Y_VALUE, 2)) ** (1 / 2)
        f = smp.sin(x) ** 3 * smp.exp(-5 * x)
        result = integrate.quad(lambda x: special.jv(2.5, x), 0, 4.5)
        print(result)


class A_I:
    def integrand(r, F_dy):
        return (-1 / np.pi) * F_dy / ((y**2 - r**2) ** (1 / 2))

    def derivative(F_y):
        return smp.diff(F_y)

    def abel_inverse_given_F_y(r, a, F_y):
        y = smp.symbols("y", real=True)
        F_dy = A_I.derivative(F_y)
        test = A_I.integrand(r, F_dy)
        print(test)
        print(smp.integrate(test, (y, r, a)))
        return smp.integrate(test, (y, r, a)).evalf()


y = smp.symbols("y", real=True)
F_y = y**3 + y
test_1 = A_I.abel_inverse_given_F_y(cv.R_VALUE, cv.A_VALUE, F_y).as_two_terms()
print(test_1)
