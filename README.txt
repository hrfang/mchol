This set of codes (gmw81, gmw1, gmw2, se90, se99, and se1) compute Cholesky
factorizations of real symmetric matrices, modified if necessary to make
them positive definite. They differ in how they compute the modification.

Given an input matrix A, each produces output matrices L, P, and E so that
            P * (A + E) * P' = L * L'
where P is a permutation matrix, L is lower triangular, and E is the
modification.

For optimization, the SE90 algorithm is not generally stable, whereas the
other 5 algorithms usually work OK.


Table of Contents
=================
- How to install?
- How to use?
- The 1+33 test matrices
- Octave/Matlab mex interface
- Modified Newton methods for unconstrained optimization


How to install?
===============
cd ./drivers
make
# it should build an archive file "./source/libmchol.a" and then executable
# "./drivers/gmw" and "./drivers/se".


How to use?
===========
Run "gmw" without arguments, and a help message should print.
It takes input/output matrices in matrix-market format (mtx file).
A Matlab script to read mtx files is here:
http://math.nist.gov/MatrixMarket/mmio/matlab/mmread.m

Examples:
---------
cd data/mtx
../../drivers/gmw -gmw81 FIRST.mtx L.mtx -P=P.mtx -E=E.mtx
# here FIRST.mtx is the input matrix.
# we should get the Cholesky factor L stored in "L.mtx",
# the permutation matrix P stored in "P.mtx", and
# the modification matrix stored in "E.mtx".


The 1+33 test matrices:
======
There are 1+33 matrices are used for testing the modified Cholesky
algorithms in the SE99 paper and our 2008 paper. The SE99 paper refers to:
------
A Revised Modified Cholesky Factorization Algorithm,
Robert B. Schnabel and Elizabeth Eskow,
SIAM J. Optim., Vol. 9 No. 4, 1135-1148, 1999.
------
The 1+33 matrices are provided in 2 format, mtx and csv, in directories
"./data/mtx/" and "./data/csv/".

The test results of the 1+33 test matrices can be found in the following
folders "./tests_mtx/" and "./tests_csv/", parallel to each other. In
"./tests_mtx/", the factorizations are performed by the executable files
under "./drivers/" and the matrix files are in mtx (matrix-market) format.
In "./tests_csv/", the factorizations are performed by the mex files under
"./mex/" and the matrix files are in csv format. See "./tests_mtx/README.txt"
and "./tests_csv/README.txt" for more information.


Octave/Matlab mex interface:
======
We provide an Octave/Matlab mex interface such that the 6 modified Cholesky
algorithms can be invoked directly in Octave/Matlab. Go to "./mex/" and
under Octave/Matlab, execute "make" for the 6 mex files of the 6 algorithms.
See "./mex/README.txt" for more information.


Modified Newton methods for unconstrained optimization:
======
In "./opt/", we provide Octave/Matlab m-files for the modified Newton
methods for unconstrained optimization with examples. It requires the mex
files under "./mex/". See "./opt/README.txt" for more information.


Questions, comments, and complaints, send us email:
Hawren Fang,        hrfang (at) yahoo (dot) com
Dianne P. O'Leary,  oleary (at) cs (dot) umd (dot) edu

Reference:
Modified Cholesky Algorithms: A Catalog with New Approaches,
Haw-ren Fang and Dianne O'Leary,
Mathematical Programming, Series A, Vol. 115, No. 2, pp. 319--349, 2008.
