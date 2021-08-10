'''
File: debug.py
Author: Yutong Dai (rothdyt@gmail.com)
File Created: 2020-08-31 21:42
Last Modified: 2021-05-28 00:05
--------------------------------------------
Description:
'''

from src.proxGradNZZ import PGSOLVER
from src.problemNZZ import Problem
from src.lossfunction import LogisticLoss, LeastSquares
from src.regularizer import GL1
from src.params import *
import src.utils as utils
import sys
import os
import numpy as np
PACKAGE_PARENT = '../'
SCRIPT_DIR = os.path.dirname(os.path.realpath(
    os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))


testType = 'logit'
if testType == 'logit':
    # datasetName = 'diabetes'
    datasetName = 'mushrooms'
    percent = 0.01
    frac = 0.1
    fileType = fileTypeDict[datasetName]
    print("Working on: {}...".format(datasetName))
    X, y = utils.set_up_xy(datasetName, fileType, dbDir='../../db')
    f = LogisticLoss(X, y, datasetName)
    p = X.shape[1]
    num_of_groups = int(frac * p)
    group = utils.gen_group(p, num_of_groups)
    Lambda = utils.lam_max(X, y, group) * percent
    r = GL1(Lambda=Lambda, group=group)
    problem = Problem(f, r)
    options = params
    options["beta"] = 1.1
    # options["printevery"] = 1
    ista = PGSOLVER(problem, "ISTA")
    result = ista.solve(X_initial=None, alpha=None,
                        tc='ITERATES', tol=1e-6, log=True, options=options)
    print("fevals:{} | gevals:{} | baks:{}".format(
        result['fevals'], result['gevals'], result['baks']))
    fista = PGSOLVER(problem, "FISTA")
    result = fista.solve(X_initial=None, alpha=None,
                         tc='ITERATES', tol=1e-6, log=True, options=options)
    print("fevals:{} | gevals:{} | baks:{}".format(
        result['fevals'], result['gevals'], result['baks']))
else:
    # ====================================================================================================
    datasetName = 'triazines_scale'
    percent = 0.1
    frac = 0.25
    fileType = fileTypeDict[datasetName]
    print("Working on: {}...".format(datasetName))
    X, y = utils.set_up_xy(datasetName, fileType, dbDir='../../db')
    f = LeastSquares(X, y, datasetName)
    p = X.shape[1]
    num_of_groups = int(frac * p)
    group = utils.gen_group(p, num_of_groups)
    Lambda = utils.lam_max(X, y, group, loss='ls') * percent
    r = GL1(Lambda=Lambda, group=group)
    problem = Problem(f, r)
    options = params
    options["beta"] = 1.1
    fista = PGSOLVER(problem, "FISTA")
    result = fista.solve(X_initial=None, alpha=None,
                         tc='ITERATES', tol=1e-6, log=True, options=options)
    print("fevals:{} | gevals:{} | baks:{}".format(
        result['fevals'], result['gevals'], result['baks']))
