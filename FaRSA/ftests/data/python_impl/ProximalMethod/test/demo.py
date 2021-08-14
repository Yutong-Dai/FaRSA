'''
File: debug.py
Author: Yutong Dai (rothdyt@gmail.com)
File Created: 2020-08-31 21:42
Last Modified: 2021-05-28 00:05
--------------------------------------------
Description:
'''
import sys
sys.path.append("../")
import numpy as np
from src.proxGradNZZ import PGSOLVER
from src.problemNZZ import Problem
from src.lossfunction import LogisticLoss, LeastSquares
from src.regularizer import GL1
from src.params import *
import src.utils as utils

params["printevery"] = 1


np.random.seed(417)
m, n = 10, 6
A = np.random.randn(m, n)
b = np.random.rand(m, 1)
for i in range(m):
    A[i, i % n:i % n+5] = 0
f = LeastSquares(A, b, 'randn106')
p = A.shape[1]
group = np.array([1, 1, 2, 2, 2, 3])
Lambda = 1.23#0.1
r = GL1(Lambda=Lambda, group=group)
problem = Problem(f, r)
options = params
options["beta"] = 1.0
ista = PGSOLVER(problem, "ISTA")
X = np.array([0.0, 0.0, 1.1, 0.0, 2.2, 3.3]).reshape(-1, 1)
result = ista.solve(X_initial=X, alpha=1.0,
                    tc='PROX', tol=1e-6, log=True, options=options)
print("fevals:{} | gevals:{} | baks:{}".format(
    result['fevals'], result['gevals'], result['baks']))
