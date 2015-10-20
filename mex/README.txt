This set of codes (gmw81, gmw1, gmw2, se90, se99, and se1) compute Cholesky
factorizations of real symmetric matrices, modified if necessary to make
them positive definite. They differ in how they compute the modification.

Given an input matrix A, each produces output matrices L, P, and E so that 
            P * (A + E) * P' = L * L'
where P is a permutation matrix, L is lower triangular, and E is the
modification.


======

In OCTAVE/MATLAB, "make" should get 6 mex files for OCTAVE/MATLAB routines
"gmw81", "gmw1", "gmw2", "se90", "se99", "se1". They are all modified
Cholesky factorizations in the form P*(A+E)*P'=L*L', where E is the
non-negative diagonal modification matrix, P is the permutation matrix for
pivoting, and L is the low triangular factor.

The usage is straightforward. for example, given a symmetric real matrix A,
------
[L,P,E] = gmw81(A);
------
gets the modified Cholesky factorization by the GMW81 algorithm. One can
also try "gmw1", "gmw2", "se90", "se99", "se1".

"first.csv" is the first matrix reported in the SE99 paper that the SE90
algorithm does not produce good result. "test_first.m" illustrates how to
apply "gmw81" to "first.csv".

The directory "../opt/" contains files for the modified Newton methods for
unconstrained optimization, requiring the mex files compiled in this
directory.
