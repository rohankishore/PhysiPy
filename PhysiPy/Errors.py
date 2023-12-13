from .constants import *
import numpy as np


class Errors:

    @staticmethod
    def error_muldiv(a, b, c, d):
        x = ((a / c) + (b / d))
        return x

    @staticmethod
    def error_addsub(a, b):
        p = (a + b)
        return p

    @staticmethod
    def percentage_error(M, E):
        print("Percentage Error:", (M / E) * 100, "%")

    @staticmethod
    def absolute_error(a, n):
        if len(a) == n:
            A = np.array(a)
            S = np.sum(A)
            E = S / n
            return E
        else:
            return ValueError

    @staticmethod
    def meanabsolute_error(a, n):
        if len(a) == n:
            b = []
            for i in a:
                b.append(abs(i))
            A = np.array(b)
            S = np.sum(A)
            M = S / n
            return M
        else:
            print("error try again")
