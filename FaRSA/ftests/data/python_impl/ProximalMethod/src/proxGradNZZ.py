'''
File: Solver.py
Author: Yutong Dai (rothdyt@gmail.com)
File Created: 2020-09-14 13:46
Last Modified: 2020-12-22 17:07
--------------------------------------------
Description:
'''
import time
import numpy as np
from src.params import params
import src.utils as utils
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.realpath(
    os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, '..')))


class PGSOLVER():
    def __init__(self, prob, method):
        self.prob = prob
        self.method = method
        self.version = "0.1 (2020-12-19)"

    def update_constant_stepsize(self, X, gradfX, alpha):
        Xtrial, zeroGroup, nonZeroGroup = self.prob.proximal(X, gradfX, alpha)
        fval_trial = self.prob.eval_fun_f(Xtrial)
        return Xtrial, fval_trial, alpha, 0, zeroGroup, nonZeroGroup

    def update_backtracking_stepsize(self, X, gradfX, fvalX, alpha, zeta, maxbak, strategy="frac", mode='once'):
        if mode == 'once':
            Xtrial, zeroGroup, nonZeroGroup = self.prob.proximal(
                X, gradfX, alpha)
            d = Xtrial - X
            bak = 0
            stepsize = 1.0
            Fold = fvalX + self.prob.eval_fun_r(X)
            while True:
                Xnew = X + stepsize * d
                fvalnew = self.prob.eval_fun_f(Xnew)
                Fnew = fvalnew + self.prob.eval_fun_r(Xnew)
                LHS = Fnew - Fold
                RHS = -1/alpha * np.sum(d * d) * stepsize * 1e-10
                sufficient_decrease = LHS <= RHS
                if sufficient_decrease:
                    if bak > 0:
                        alpha *= 0.8
                    return Xnew, fvalnew, alpha, bak, zeroGroup, nonZeroGroup
                stepsize *= 0.5
                bak += 1
        else:
            bak = 0
            while True:
                Xtrial, zeroGroup, nonZeroGroup = self.prob.proximal(
                    X, gradfX, alpha)
                fval_trial = self.prob.eval_fun_f(Xtrial)
                d = Xtrial - X
                dirder = np.dot(gradfX.T, d)[0][0]
                d_sq = np.dot(d.T, d)[0][0]
                LHS = fval_trial - fvalX
                RHS = dirder + 0.5 * (1 / alpha) * d_sq
                if LHS <= RHS:
                    return Xtrial, fval_trial, alpha, bak, zeroGroup, nonZeroGroup
                if bak == maxbak:
                    return None, None, None, bak
                if strategy == "frac":
                    alpha *= zeta
                elif strategy == "model":
                    L_local = 2 * (LHS - dirder) / d_sq
                    alpha = max(alpha * zeta, 1 / L_local)
                bak += 1

    def set_init_alpha(self, x):
        s = 1e-2
        _ = self.prob.eval_fun_f(x)
        gradfx = self.prob.eval_grad_f(x)
        y = x - s * gradfx
        while True:
            if utils.l2_norm(y - x) > 1e-8:
                _ = self.prob.eval_fun_f(y)
                gradfy = self.prob.eval_grad_f(y)
                alpha = utils.l2_norm(
                    x - y) / (1 * utils.l2_norm(gradfx - gradfy))
                break
            else:
                s *= 10
                y = x - s * gradfx
        return alpha

    def check_termination(self, tc, X, Xtrial, tol, bak, maxbak, time_so_far, maxtime, iteration, maxiters):
        if tc == 'ITERATES':
            diff = Xtrial - X
            optim = np.sqrt(np.dot(diff.T, diff))[0][0]
            X_trial_norm = np.sqrt(np.dot(Xtrial.T, Xtrial))[0][0]
            optim = optim / X_trial_norm
            termination = (optim <= tol * max(1, 1 / X_trial_norm))
            code = 0
        elif tc == "PROX":
            diff = Xtrial - X
            optim = np.sqrt(np.dot(diff.T, diff))[0][0]
            termination = optim <= tol * max(self.s0_norm, 1)
            code = 0
        else:
            raise ValueError("Invalid termination conditions!")
        if time_so_far >= maxtime:
            termination = True
            code = 2
        if iteration == maxiters:
            termination = True
            code = 1
        if bak >= maxbak:
            termination = True
            code = -1
        return optim, termination, code

    def get_group_structure(self, X):
        nz = 0
        for i in range(self.prob.K):
            start, end = self.prob.starts[i], self.prob.ends[i]
            X_Gi = X[start:end]
            if (np.sum(np.abs(X_Gi)) == 0):
                nz += 1
        nnz = self.prob.K - nz
        return nnz, nz

    def solve(self, X_initial=None, alpha=None, tc='PROX', tol=1e-6, log=True, options=None):
        print("Solver:{} | Version:{}".format(self.method, self.version))
        if not options:
            options = params
            print("Default options:")
        zeta, maxbak, maxtime, maxiters = options['zeta'], options[
            'maxback'], options['max_time'], options['max_iter']
        stepsizeRule = options['stepsizeRule']
        stepsizeStrategy = options['stepsizeStrategy']
        print(" tol:{} | maxtime:{} | maxiter:{}".format(tol, maxtime, maxiters))
        print(" beta:{:3.3e} | zeta:{:3.3e} | eta:{:3.3e}".format(
            options['beta'], zeta, options['eta']))
        print(" stepsizeRule:{} | stepsizeStrategy:{}".format(
            stepsizeRule, stepsizeStrategy))
        try:
            outID = self.prob.f.datasetName
        except AttributeError:
            outID = None
        start_time = time.time()
        if X_initial is None:
            X = np.zeros([self.prob.p, 1])
            provide_initial = False
        else:
            X = X_initial
            provide_initial = True
        if alpha is None:
            print(' alpha is [None]; Set up by {}'.format(
                "LipschitzEstimation"))
            alpha = self.set_init_alpha(X)
        else:
            print(' proxStepsize is [{:3.3e}]'.format(alpha))
        result = {}
        iteration = 0
        fevals = 0
        gevals = 0
        baks = 0
        # collect information at the starting point
        fvalX = self.prob.eval_fun_f(X)
        gradfX = self.prob.eval_grad_f(X)
        FvalX = fvalX + self.prob.eval_fun_r(X)
        optim = "-"
        nnz, nz = self.get_group_structure(X)
        fevals += 1
        gevals += 1
        # extra setups required by the FISTA
        if self.method == "FISTA":
            y = X
            fvaly = fvalX
            gradfy = gradfX
            t = 1
            alpha_old = alpha
            tfocs_t = options['tfocs_t']
            print(" use tfocs t seq:{}".format(tfocs_t))
        extra_time = 0
        if log:
            utils.print_problem(self.method, provide_initial,
                                self.version, self.prob.r.Lambda, self.prob.K, outID)
            utils.print_header_PG_NZZ(outID)
        while True:
            if log and iteration % options['printevery'] == 0:
                utils.print_iterates_PG_NZZ(
                    iteration, FvalX, optim, nnz, nz, outID)
            if stepsizeRule == "bkt":
                if self.method == "ISTA":
                    Xtrial, fvalXtrial, alpha, bak, zeroGroup, nonZeroGroup = self.update_backtracking_stepsize(X, gradfX, fvalX,
                                                                                                                alpha, zeta, maxbak, stepsizeStrategy)
                elif self.method == "FISTA":
                    Xtrial, fvalXtrial, alpha, bak, zeroGroup, nonZeroGroup = self.update_backtracking_stepsize(y, gradfy, fvaly,
                                                                                                                alpha, zeta, maxbak, stepsizeStrategy)
                else:
                    raise ValueError("Invalid method!")
                # for "PROX" termination purpose
                if iteration == 0:
                    s0 = Xtrial - X
                    self.s0_norm = np.sqrt(np.dot(s0.T, s0))[0][0]
            elif stepsizeRule == "cst":
                if self.method == "ISTA":
                    Xtrial, fvalXtrial, alpha, bak, zeroGroup, nonZeroGroup = self.update_constant_stepsize(
                        X, gradfX, alpha)
                elif self.method == "FISTA":
                    Xtrial, fvalXtrial, alpha, bak, zeroGroup, nonZeroGroup = self.update_constant_stepsize(
                        y, gradfy, alpha)
                else:
                    raise ValueError("Invalid method!")
            else:
                raise ValueError("Invalid stepsizeRule!")
            if log and iteration % options['printevery'] == 0:
                utils.print_update_PG(bak, alpha, outID)
            iteration += 1
            FvalXtrial = fvalXtrial + self.prob.eval_fun_r(Xtrial)
            nz, nnz = len(zeroGroup), len(nonZeroGroup)
            # (# feval in linesearch = # backtrack + 1)
            fevals += bak + 1
            baks += bak
            time_so_far = time.time() - start_time
            optim, termination, code = self.check_termination(
                tc, X, Xtrial, tol, bak, maxbak, time_so_far, maxtime, iteration, maxiters)
            # terminate algorithm and collect data
            if termination:
                total_time = time.time() - start_time
                result = {
                    'X': Xtrial, 'iteration': iteration, 'time': total_time, 'F': FvalXtrial,
                    'nz': nz, 'nnz': nnz, 'status': code,
                    'fevals': fevals, 'gevals': gevals, 'baks': baks, 'extra_time': extra_time
                }
                if log:
                    if iteration % options['printevery'] != 0:
                        utils.print_iterates_PG_NZZ(
                            iteration, FvalXtrial, optim, self.prob.K - nz, nz, outID)
                    utils.print_exit(code, outID)
                    utils.print_result_NZZ(result, outID)
                print("Time Elapsed:{:5.3f} | Iterations:{} | Extrapolation time:{:5.3f}\n".format(
                    total_time, iteration, extra_time))
                print("-" * 50)
                return result
            # prepare for the next iteration
            if self.method == "FISTA":
                start = time.time()
                if tfocs_t:
                    ratio = alpha_old / alpha
                    t_new = (1 + np.sqrt(1 + 4 * ratio * t**2)) / 2
                    alpha_old = alpha
                else:
                    t_new = (1 + np.sqrt(1 + 4 * t**2)) / 2
                y = Xtrial + ((t - 1) / t_new) * (Xtrial - X)
                t = t_new
                fvaly = self.prob.eval_fun_f(y)
                fevals += 1
                extra_time += time.time() - start
                gradfy = self.prob.eval_grad_f(y)
                gevals += 1
                X = Xtrial
                FvalX = FvalXtrial
            else:
                X = Xtrial
                fvalX = fvalXtrial
                FvalX = FvalXtrial
                gradfX = self.prob.eval_grad_f(X)
                gevals += 1
            alpha *= options['beta']
