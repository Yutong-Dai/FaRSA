FaRSA
=====

This software is under development.  Please check back later.

Overview
--------

FaRSA (Fast Reduced Space Algorithm) is a software package for solving minimization problems.  It is designed to locate a minimizer of

```
min     f(x) + r(x)
x ∈ Rⁿ
```

where ```f : Rⁿ -> R``` is continuously differentiable and ```r : Rⁿ -> R``` is convex.  In particular, it is intended for problems in which ```f``` is a data-fitting function and ```r``` is a (group) sparsity-inducing regularization function.

FaRSA is written in C++ and is released under the ??? License.  The main authors are [Frank E. Curtis](http://coral.ise.lehigh.edu/frankecurtis/) and [Daniel P. Robinson](https://coral.ise.lehigh.edu/danielprobinson/).  For a list of all contributors, please see the [AUTHORS file](FaRSA/AUTHORS).

Compiling FaRSA requires BLAS and LAPACK.  The code for these packages is not provided in this repository, which are available under different conditions and licenses than those for FaRSA.

Citing FaRSA
------------

FaRSA is provided free of charge so that it might be useful to others.  Please send e-mail to [Frank E. Curtis](http://coral.ise.lehigh.edu/frankecurtis/) and [Daniel P. Robinson](https://coral.ise.lehigh.edu/danielprobinson/) with success stories or other feedback.  If you use FaRSA in your research, then please cite the following papers:

- F. E. Curtis, Y. Dai, and D. P. Robinson. "A Subspace Acceleration Method for Minimization Involving a Group Sparsity-Inducing Regularizer." [arXiv:2007.14951](https://arxiv.org/abs/2007.14951), 2020.
- T. Chen, F. E. Curtis, and D. P. Robinson. "FaRSA for l1-Regularized Convex Optimization: Local Convergence and Numerical Experience." Optimization Methods and Software, 33(2):396–415, 2018.
- T. Chen, F. E. Curtis, and D. P. Robinson. "A Reduced-Space Algorithm for Minimizing $\ell_1$-Regularized Convex Functions." SIAM Journal on Optimization, 27(3):1583–1610, 2017.
