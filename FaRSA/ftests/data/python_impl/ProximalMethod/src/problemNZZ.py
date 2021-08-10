'''
File: problem.py
Author: Yutong Dai (rothdyt@gmail.com)
File Created: 2020-12-19 22:46
Last Modified: 2020-12-22 12:12
--------------------------------------------
Description:
'''
from numba import jit
import numpy as np


# @jit(nopython=True, cache=True)
def _proximal_gradient_jit(X, alpha, gradf, K, p, starts, ends, Lambda_group):
    proximal = np.zeros((p, 1))
    nonZeroGroup = []
    zeroGroup = []
    for i in range(K):
        start, end = starts[i], ends[i]
        XG_i = X[start:end]
        gradfG_i = gradf[start:end]
        gradient_step = XG_i - alpha * gradfG_i
        gradient_step_norm = np.sqrt(np.dot(gradient_step.T, gradient_step))[0][0]
        if gradient_step_norm != 0:
            temp = 1 - ((Lambda_group[i] * alpha) / gradient_step_norm)
        else:
            temp = -1
        if temp > 0:
            nonZeroGroup.append(i)
        else:
            zeroGroup.append(i)
        proximal[start:end] = max(temp, 0) * gradient_step
    return proximal, zeroGroup, nonZeroGroup


class Problem:
    def __init__(self, f, r):
        self.f = f
        self.r = r
        self.K = self.r.K
        self.n, self.p = self.f.n, self.f.p
        self.starts = self.r.starts
        self.ends = self.r.ends
        self.Lambda_group = self.r.Lambda_group

    def eval_fun_f(self, X):
        return self.f.evaluate_function_value(X)

    def eval_grad_f(self, X=None):
        return self.f.gradient()

    def eval_fun_r(self, X):
        return self.r.eval_f_vec(X)

    def proximal(self, X, gradfX, alpha):
        Xnew, zeroGroup, nonZeroGroup = _proximal_gradient_jit(X, alpha, gradfX,
                                                               self.K, self.p, self.starts,
                                                               self.ends, self.Lambda_group)
        return Xnew, zeroGroup, nonZeroGroup
