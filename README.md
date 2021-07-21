# A C++ Implementation of Modified Cholesky Factorizations

This set of codes compute Cholesky factorizations of real symmetric matrices,
modified if necessary to make them positive definite.
We have included implementations of modified Cholesky algorithms
**gmw81**, **gmw1**, **gmw2**, **se90**, **se99**, and **se1**.
They differ in how they compute the modification.

To be precise, given a symmetric matrix <span class="math display"><em>A</em></span>, each algorithm produces
output matrices <span class="math display"><em>L</em></span>,
<span class="math display"><em>P</em></span>, and
<span class="math display"><em>E</em></span>, so that

<span class="math display"><em>P</em>(<em>A</em>+<em>E</em>)<em>P</em><sup><em>T</em></sup> = <em>L</em><em>L</em><sup><em>T</em></sup>,</span>

where \\(P\\) is a permutation matrix for pivoting, $$L$$ is lower triangular,
and \\(E\\) is the modification.

In practice, for non-convex optimization, the **se90** algorithm is not generally stable,
whereas the other 5 algorithms usually work fine.


## Table of Contents

- [Why modified Cholesky factorizations?](#why-modified-cholesky-factorizations)
- [How to install?](#how-to-install)
- [How to use?](#how-to-use)
- [The 1+33 test matrices](#the-133-test-matrices)
- [Matlab/Octave mex interface](#matlaboctave-mex-interface)
- [Modified Newton methods for unconstrained optimization](#modified-newton-methods-for-unconstrained-optimization)


## Why modified Cholesky factorizations?

To minimize a function $f:\mathbb{R}^n\rightarrow\mathbb{R}$ by Newton's method,
we solve a linear system $Ax=b$ for the search direction at every iteration,
where $A\in\mathbb{R}^{n\times n}$ is the Hessian matrix and $b\in\mathbb{R}^n$ is the negated gradient.

The Hessian matrix $A$ is symmetric
under the assumption that second partial derivatives
of the objective function $f$ are continuous.
Moreover, $A$ is symmetric positive definite (SPD) and the search
direction from solving $Ax=b$ is descent,
if the objective function $f$ is strictly convex.
Coupled with a globalization technique such as a line search or trust region
technique, Newton's method is guaranteed to converge to the minimum,
which is unique and hence global due to convexity.

When the objective function $f$ is not convex, the Hessian $A$ may
be indefinite, and Netwon's direction may not be descent.
As a result, convergence to a minimum is no longer guaranteed.
A modified Cholesky algorithm perturbs $A$ to be a positive definite
$A+E$ to ensure a descent search direction and therefore the convergence
to a local minimum.


## How to install?

```
$ cd ./drivers
$ make
```

It should build an archive file `./source/libmchol.a` and then executable
`./drivers/gmw` and `./drivers/se`.


## How to use?

Run `gmw` without arguments, and a help message should print.
It takes input/output matrices in matrix-market format (mtx file).

Here is an example.
```
$ cd data/mtx
$ ../../drivers/gmw -gmw81 FIRST.mtx L.mtx -P=P.mtx -E=E.mtx
```

Here `FIRST.mtx` is the input matrix.
We should get the Cholesky factor $L$ stored in `L.mtx`,
the permutation matrix $P$ stored in `P.mtx`, and
the modification matrix stored in `E.mtx`.

For Matlab/Octave users, a Matlab/Octave script to read mtx files is here:

[http://math.nist.gov/MatrixMarket/mmio/matlab/mmread.m][mmread_link]

[mmread_link]: http://math.nist.gov/MatrixMarket/mmio/matlab/mmread.m



## The 1+33 test matrices

There are 1+33 matrices are used for testing the modified Cholesky
algorithms in the **se99** paper and
our 2008 paper (see [reference](#reference)).
The **se99** paper refers to:

Robert B. Schnabel and Elizabeth Eskow,
"A Revised Modified Cholesky Factorization Algorithm,"
SIAM J. Optim., Vol. 9 No. 4, 1135-1148, 1999.

The 1+33 matrices are provided in 2 format, mtx and csv, in directories
`./data/mtx/` and `./data/csv/`.

The test results of the 1+33 test matrices can be found in the following
folders `./tests_mtx/` and `./tests_csv/`, parallel to each other. In
`./tests_mtx/`, the factorizations are performed by the executable files
under `./drivers/` and the matrix files are in mtx (matrix-market) format.
In `./tests_csv/`, the factorizations are performed by the mex files under
`./mex/` and the matrix files are in csv format. See `./tests_mtx/README.txt`
and `./tests_csv/README.txt` for more information.


## Matlab/Octave mex interface

We provide an Matlab/Octave mex interface such that the 6 modified Cholesky
algorithms can be invoked directly in Matlab/Octave. Go to `./mex/` and
under Matlab/Octave, execute `make` for the 6 mex files of the 6 algorithms.
See `./mex/README.txt` for more information.


## Modified Newton methods for unconstrained optimization

In `./opt/`, we provide Matlab/Octave m-files for the modified Newton
methods for unconstrained optimization with examples. It requires the mex
files under `./mex/`. See `./opt/README.txt` for more information.


## Contact us

Questions, comments, and complaints, send us email:
- Hawren Fang,        hrfang (at) yahoo (dot) com
- Dianne P. O'Leary,  oleary (at) cs (dot) umd (dot) edu


## Reference

Haw-ren Fang and Dianne O'Leary,
"Modified Cholesky Algorithms: A Catalog with New Approaches,"
Mathematical Programming, Series A, Vol. 115, No. 2, pp. 319--349, 2008.
