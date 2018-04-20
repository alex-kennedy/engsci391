import numpy as np

class LinearProgram(object):
    def __init__(cost, a_matrix, b):
        self.c = c
        self.a = a
        self.b = rhs

        mn = a_matrix.shape
        self.m = mn[0]
        self.n = mn[0]

        self.basic_vars = np.arange(self.m + 1, self.m + self.n)
        self.b_inv = np.identity(self.m)
        self.phase1 = True

    def get_pi(self):
        return


    def find_ev(self):
        return


    def find_rm(self):
        return

    def update(self):
        return
