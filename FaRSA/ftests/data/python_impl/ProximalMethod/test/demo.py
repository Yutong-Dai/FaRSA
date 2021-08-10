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

# testType = 'logit'
# if testType == 'logit':
#     # datasetName = 'diabetes'
#     datasetName = 'mushrooms'
#     percent = 0.01
#     frac = 0.1
#     fileType = fileTypeDict[datasetName]
#     print("Working on: {}...".format(datasetName))
#     X, y = utils.set_up_xy(datasetName, fileType,
#                            dbDir='/Users/ym/Documents/GroupFaRSA/db')
#     f = LogisticLoss(X, y, datasetName)
#     p = X.shape[1]
#     num_of_groups = int(frac * p)
#     group = utils.gen_group(p, num_of_groups)
#     Lambda = utils.lam_max(X, y, group) * percent
#     r = GL1(Lambda=Lambda, group=group)
#     problem = Problem(f, r)
#     options = params
#     options["beta"] = 1.1
#     # options["printevery"] = 1
#     ista = PGSOLVER(problem, "ISTA")
#     result = ista.solve(X_initial=None, alpha=None,
#                         tc='ITERATES', tol=1e-6, log=True, options=options)
#     print("fevals:{} | gevals:{} | baks:{}".format(
#         result['fevals'], result['gevals'], result['baks']))
#     fista = PGSOLVER(problem, "FISTA")
#     result = fista.solve(X_initial=None, alpha=None,
#                          tc='ITERATES', tol=1e-6, log=True, options=options)
#     print("fevals:{} | gevals:{} | baks:{}".format(
#         result['fevals'], result['gevals'], result['baks']))
# else:
#     # ====================================================================================================
#     datasetName = 'triazines_scale'
#     percent = 0.1
#     frac = 0.25
#     fileType = fileTypeDict[datasetName]
#     print("Working on: {}...".format(datasetName))
#     X, y = utils.set_up_xy(datasetName, fileType, dbDir='../../db')
#     f = LeastSquares(X, y, datasetName)
#     p = X.shape[1]
#     num_of_groups = int(frac * p)
#     group = utils.gen_group(p, num_of_groups)
#     Lambda = utils.lam_max(X, y, group, loss='ls') * percent
#     r = GL1(Lambda=Lambda, group=group)
#     problem = Problem(f, r)
#     options = params
#     options["beta"] = 1.1
#     fista = PGSOLVER(problem, "FISTA")
#     result = fista.solve(X_initial=None, alpha=None,
#                          tc='ITERATES', tol=1e-6, log=True, options=options)
#     print("fevals:{} | gevals:{} | baks:{}".format(
#         result['fevals'], result['gevals'], result['baks']))
