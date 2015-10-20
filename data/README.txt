This set of codes (gmw81, gmw1, gmw2, se90, se99, and se1) compute Cholesky
factorizations of real symmetric matrices, modified if necessary to make
them positive definite. They differ in how they compute the modification.

Given an input matrix A, each produces output matrices L, P, and E so that 
            P * (A + E) * P' = L * L'
where P is a permutation matrix, L is lower triangular, and E is the
modification.


======

The folder 33matrices contains the 33 matrices from Betty Eskow that were
used in the experiments of the SE99 paper. The first matrix identified for
which SE90 has a problem is also included. There are 34 matrices in total.

Two formats are provided: mtx (matrix-market) and csv, in directories
"./mtx/" and "./csv/", respectively.
