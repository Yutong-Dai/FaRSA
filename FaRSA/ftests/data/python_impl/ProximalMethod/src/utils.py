'''
File: utils.py
Author: Yutong Dai (rothdyt@gmail.com)
File Created: 2020-09-14 12:32
Last Modified: 2020-12-22 17:10
--------------------------------------------
Description:
'''
import numpy as np
from sklearn.datasets import load_svmlight_file


def set_up_xy(datasetName, fileType='txt', dbDir='../db', to_dense=False, order=None):
    filepath = dbDir + "/{}.{}".format(datasetName, fileType)
    data = load_svmlight_file(filepath)
    X, y = data[0], data[1].reshape(-1, 1)
    # if datasetName in ['gisette']:
    #     to_dense = True
    if to_dense:
        print("  Begin converting {}...".format(datasetName))
        X = X.toarray(order)
        print("  Finish converting!")
    return X, y


def l2_norm(x):
    return np.sqrt(np.dot(x.T, x))[0][0]


# def lam_max(X, y, group):
#     """
#     Reference: Yi Yang and Hui Zou. A fast unified algorithm for solving group-lasso penalize learning problems. Page 22.
#     """
#     # beta_1 = np.log(np.sum(y == 1) / np.sum(y == -1))
#     beta_1 = 0
#     ys = y / (1 + np.exp(beta_1 * y))
#     nabla_L = X.T@ys / X.shape[0]
#     # print(beta_1, norm(ys, 2)**2, nabla_L.T)
#     lam_max = -1
#     unique_groups, group_frequency = np.unique(group, return_counts=True)
#     K = len(unique_groups)
#     for i in range(K):
#         sub_grp = nabla_L[group == (i + 1)]
#         temp = l2_norm(sub_grp) / np.sqrt(group_frequency[i])
#         if temp > lam_max:
#             lam_max = temp
#     return lam_max
def lam_max(X, y, group, loss='logit'):
    """
    Reference: Yi Yang and Hui Zou. A fast unified algorithm for solving group-lasso penalize learning problems. Page 22.
    """
    # beta_1 = np.log(np.sum(y == 1) / np.sum(y == -1))
    beta_1 = np.zeros((X.shape[1], 1))
    lam_max = -1
    unique_groups, group_frequency = np.unique(group, return_counts=True)
    K = len(unique_groups)
    if loss == 'logit':
        ys = y / (1 + np.exp(y * (X @ beta_1)))
        nabla_L = X.T @ ys / X.shape[0]
    elif loss == 'ls':
        ys = y - X @ beta_1
        nabla_L = (ys.T @ X).T / X.shape[0]
    else:
        raise ValueError("Invalid loss!")
    for i in range(K):
        sub_grp = nabla_L[group == (i + 1)]
        temp = l2_norm(sub_grp) / np.sqrt(group_frequency[i])
        if temp > lam_max:
            lam_max = temp
    return lam_max


def gen_group(p, K):
    dtype = type(K)
    if dtype == int:
        group = K * np.ones(p)
        size = int(np.floor(p / K))
        for i in range(K):
            start_ = i * size
            end_ = start_ + size
            group[start_:end_] = i + 1
    elif dtype == np.ndarray:
        portion = K
        group = np.ones(p)
        chunk_size = p * portion
        start_ = 0
        for i in range(len(portion)):
            end_ = start_ + int(chunk_size[i])
            group[start_:end_] = i + 1
            start_ = int(chunk_size[i]) + start_
    return group


def print_problem(method, provide_initial, version, Lambda, K, outID=None):
    if outID is not None:
        filename = '{}.txt'.format(outID)
    else:
        filename = 'log.txt'
    with open(filename, "a") as logfile:
        contents = "\n" + "=" * 80
        contents += "\n          Proximal Gradient Type Method   (version:{})  \n".format(version)
        contents += "=" * 80
        contents += "\nMethod:{:.>48}{}\n".format("", method)
        contents += "Provide Initial:{:.>39}{}\n".format("", provide_initial)
        contents += "Penalty Parameter:{:.>37}{:3.4e}\n".format('', Lambda)
        contents += "Number of groups:{:.>40}\n".format(K)
        contents += "\n"
        logfile.write(contents)


def print_header_PG(outID=None):
    if outID is not None:
        filename = '{}.txt'.format(outID)
    else:
        filename = 'log.txt'
    column_titles = '  Iter            F            optim        bak     alpha\n'
    with open(filename, "a") as logfile:
        logfile.write(column_titles)


def print_iterates_PG(iteration, F, optim, outID=None):
    if outID is not None:
        filename = '{}.txt'.format(outID)
    else:
        filename = 'log.txt'
    if isinstance(optim, str):
        contents = " {:5d}       {:3.6e}   ----------- |".format(iteration, F)
    else:
        contents = " {:5d}       {:3.6e}   {:3.5e} |".format(iteration, F, optim)
    with open(filename, "a") as logfile:
        logfile.write(contents)


def print_update_PG(bak, alpha, outID=None):
    if outID is not None:
        filename = '{}.txt'.format(outID)
    else:
        filename = 'log.txt'
    contents = " {:3d}    {:3.4e}   |\n".format(bak, alpha)
    with open(filename, "a") as logfile:
        logfile.write(contents)


def print_exit(status, outID=None):
    if outID is not None:
        filename = '{}.txt'.format(outID)
    else:
        filename = 'log.txt'
    contents = '\n' + "=" * 30 + '\n'
    if status == -1:
        contents += 'Exit: Line Search Failed\n'
    elif status == 0:
        contents += 'Exit: Optimal Solution Found\n'
    elif status == 1:
        contents += 'Exit: Iteration limit reached\n'
    elif status == 2:
        contents += 'Exit: Time limit reached\n'
    with open(filename, "a") as logfile:
        logfile.write(contents)
        print(contents)


def print_result(results, outID):
    if outID is not None:
        filename = '{}.txt'.format(outID)
    else:
        filename = 'log.txt'
    contents = "\nSummary\n"
    contents += "=" * 50
    contents += "\nFinal Objective Value F:{:.>28}{:3.6e}\n".format("", results['Fval'])
    contents += "Total Iterations:{:.>35}{}\n".format("", results['iters'])
    contents += "Total time:{:.>41}{:5.3e}\n\n".format("", results['time'])
    with open(filename, "a") as logfile:
        logfile.write(contents)


def print_header_PG_NZZ(outID=None):
    if outID is not None:
        filename = '{}.txt'.format(outID)
    else:
        filename = 'log.txt'
    column_titles = '  Iter            F            optim        nnz        nz    bak     alpha\n'
    with open(filename, "a") as logfile:
        logfile.write(column_titles)


def print_iterates_PG_NZZ(iteration, F, optim, nz, z, outID=None):
    if outID is not None:
        filename = '{}.txt'.format(outID)
    else:
        filename = 'log.txt'
    if isinstance(optim, str):
        contents = " {:5d}       {:3.6e}   -----------   {:5d}     {:5d} |".format(iteration, F, nz, z)
    else:
        contents = " {:5d}       {:3.6e}   {:3.5e}   {:5d}     {:5d} |".format(iteration, F, optim, nz, z)
    with open(filename, "a") as logfile:
        logfile.write(contents)


def print_result_NZZ(results, outID):
    if outID is not None:
        filename = '{}.txt'.format(outID)
    else:
        filename = 'log.txt'
    contents = "\nSummary\n"
    contents += "=" * 80
    contents += "\nFinal Objective Value F:{:.>28}{:3.6e}\n".format("", results['F'])
    contents += "Total Iterations:{:.>35}{}\n".format("", results['iteration'])
    contents += "Total time:{:.>41}{:5.3e}\n".format("", results['time'])
    contents += "#Zero Groups:{:.>39}{}\n".format("", results['nz'])
    contents += "Extrapolation time:{:.>33}{:5.3e}\n".format("", results['extra_time'])
    contents += "Function evaluations:{:.>31}{:d}\n".format("", results['fevals'])
    contents += "Gradient evaluations:{:.>31}{:d}\n".format("", results['gevals'])
    contents += "Backtrackings:{:.>38}{:d}\n\n".format("", results['baks'])
    with open(filename, "a") as logfile:
        logfile.write(contents)
