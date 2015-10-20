This set of codes (gmw81, gmw1, gmw2, se90, se99, and se1) compute Cholesky
factorizations of real symmetric matrices, modified if necessary to make
them positive definite. They differ in how they compute the modification.

Given an input matrix A, each produces output matrices L, P, and E so that 
            P * (A + E) * P' = L * L'
where P is a permutation matrix, L is lower triangular, and E is the
modification.


======

This directory contains the source files of MCHOL.

On Linux/Unix systems, "make" should produce the achieve file "libmchol.a".
In case of errors, try tweaking "Makefile".

For binary executable drivers, go to the directory "../drivers" and see the
README.txt therein.
