# 99problems

Do `make install` to setup Python virtual environment. Activate the virtual environment with `source venv/bin/activate`. Do `pytest test` to run verification examples. Further documentation to follow.

Adapted and developed from work done with [artofscience](https://github.com/artofscience); in particular, linear programming and adaptive move-limits on self-weight topology optimisation.

The GCMMA is taken from [arjendeetman](https://github.com/arjendeetman/GCMMA-MMA-Python), November 2022.

From [SOAPs](https://github.com/dirkmunro89/SOAPs):

## Description

This repository contains a simple procedural representation of a general optimization (nonlinear programming) algorithm. In particular, it represents the established form of nonlinear programming algorithm found in the field of structural optimization, the particularities of which is in development since at least the 1970´s [[2]](#1). That being said, it is not incorrect to draw the equivalence with Newton's method (or the Newton-Raphson method). Apparently, as early as 1740, Thomas Simpson (of numerical integration fame) described Newton's method as an iterative (as in incremental) method for solving systems of nonlinear equations, noting that it can be applied to optimization problems by solving for the stationary point [[3]](#2). Today, the form of nonlinear programming algorithm found in structural optimization is referred to *Sequential Approximate Optimization (SAO)*, and it is typically characterized by:

- *positive* (and bounded) *design variables*, representing physical quantities such as *e.g.* structural member thickness, cross-sectional area, or material density;

- *inequality constraints*, representing *e.g.* design restrictions, a restriction on overall material usage (a material distribution problem), or such quantities as the maximum displacement, or stresses permitted in the structure;

- repeated evaluation of computationally expensive *structural analysis* routines (often external finite element modeling packages); that is *simulation-based*, in general, insofar as *simulation* refers to the modeling of a physical system via numerical solution of a set of partial differential equations (PDEs);

- utilization of *intervening variables (variable transformations) and function approximation concepts* to arrive at computationally tractable yet reasonably convergent approximate subproblems;

- and *dual subproblem formulations* (e.g. Falk [[4]](#3)), permitting simple numerical solution by readily available bounded function minimization routines (e.g. L-BFGS-B). (This is true in a traditional sense, at least. Primal-dual (*e.g.* interior-point) methods are increasingly available, at the expense of the elegance and the simplicity of the classic dual formulations.)

It should however be noted that the nonlinear programming paradigm represented here, is not in general restricted to: *(i)* inequality constraints [[5]](#4), *(ii)* positive (and natural) design variables (see *e.g.* simultaneous analysis and design (SAND) [[6]](#5) and Space-Time [[7]](#6) formulations), and *(iii)* as already mentioned, dual subproblem formulations [[8]](#7).

## Keywords (in work)

Sequential Approximate Optimization (SAO)
Method of Moving Asymptotes (MMA)
CONLIN
Sequential Approximate Optimization based on Incomplete Taylor Series expansions (SAOi)
Quadratically Constrained Quadratic Programming
Quadratic Programming (linear constraints)
Dual method (Falk)
Bayesian Optimization Global

## References
<a id="9">[1]</a>
Svanberg, K. (1987)
The method of moving asymptotes&mdash;a new method for structural optimization.
International Journal for Numerical Methods in Engineering, 24:359–373.

<a id="1">[2]</a>
Fleury, C. (1979)
Structural weight optimization by dual methods of convex programming.
International Journal for Numerical Methods in Engineering, 14(12):1761–1783.

<a id="2">[3]</a>
Wikipedia (2002).
[Newton's method](https://en.wikipedia.org/wiki/Newton%27s_method).

<a id="3">[4]</a>
Falk, J.E. (1967)
Lagrange multipliers and nonlinear programming.
Journal of Mathematical Analysis and Applications, 19(1):141–159.

<a id="4">[5]</a>
Cronje, M. *et al.* (2019)
Brief note on equality constraints in pure dual SAO settings.
Structural and Multidisciplinary Optimization, 59(5):1853–1861.

<a id="5">[6]</a>
Haftka, R.T. (1985) 
Simultaneous analysis and design. 
American Institute of Aeronautics and Astronautics Journal, 23(7):1099–1103.

<a id="6">[7]</a>
Wu, J. (2020) 
[Space-time topology optimization for additive manufacturing](https://doi.org/10.1007/s00158-019-02420-6). 
Structural and Multidisciplinary Optimization, 61:1-18.

<a id="7">[8]</a>
Zillober, C. (2001) 
A combined convex approximation&mdash;interior point approach for large scale
nonlinear programming. 
Optimization and Engineering, 2(1):51–73.
