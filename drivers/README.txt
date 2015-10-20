This set of codes (gmw81, gmw1, gmw2, se90, se99, and se1) compute Cholesky
factorizations of real symmetric matrices, modified if necessary to make
them positive definite. They differ in how they compute the modification.

Given an input matrix A, each produces output matrices L, P, and E so that 
            P * (A + E) * P' = L * L'
where P is a permutation matrix, L is lower triangular, and E is the
modification.

======
We need the achieve file "../sources/libmchol.a". If it is not yet
available, go to "../sources" and "make" for it. See also
"../sources/README.txt".

After obtaining "../sources/libmschol.a". "make" in this directory should
get two binary executable files "gmw" and "se".

Run "gmw" with arguments and get the help message:
======
Usage: ./gmw INPUT_MATRIX_A OUTPUT_MATRIX_L [OPTION] [OPTION] ...
  Modified Cholesky algorithms in the GMW family.
  Given the input symmtric matrix A, the factorization is in the form
                P*(A+E)*P^T = L*L^T,
  where L is the lower triangular modified Cholesky factor, E is the
  modification matrix which is diagonal and nonnegative, and P is the
  permutation matrix from pivoting.
  The input/output matrices are dense and stored in Matrix-Market format.
Arguments:
  INPUT_MATRIX_A is the file name of the input symmetric matrix A.
  OUTPUT_MATRIX_L is the file name of the output modified Cholesky factor L.
  The log will print out the permutation corresponding to P, as well as the
  the diagonal of P*E*P^T (the modification).
Optional parameters:
  -P=MATRIX_FILE      save the permutation matrix P to file MATRIX_FILE
  -E=MATRIX_FILE      save the diagonal modification matrix E to MATRIX_FILE
  -help               display this help and exit
  -gmw81              for the GMW81  algorithm
  -gmw1               for the GMW-I  algorithm
  -gmw2               for the GMW-II algorithm
                      if "-gmw81",  "-gmw1", or "-gmw2" is specified,
                      then the following parameter specification is ignored
  -delta=VALUE        the modification tolerance delta
  -pivot=INDEX        the pivoting method (default is no pivoting)
                      1 for pivoting by maximum diagonal element,
                      2 for pivoting by maximum diagonal magnitude
  -type2              for a type II algorithm (default is type I)
  -nondecreaseing     for the diagonal modifications being non-decreasing
                      (the default allows decreasing modifications)
  -relax=VALUE        relaxation factor, which must be positive and <= 1.0
                      (the default is not to use the relaxation strategy)
  -special_last       for invoking the special treatment for the last 1-by-1
                      or 2-by-2 steps used in the SE99 algorithm
                      (the default is not to use this special treatment)


Run "se" with arguments and get the help message:
======
Usage: ./se INPUT_MATRIX_A OUTPUT_MATRIX_L [OPTION] [OPTION] ...
  Modified Cholesky algorithms in the SE family.
  Given the input symmtric matrix A, the factorization is in the form
                P*(A+E)*P^T = L*L^T,
  where L is the lower triangular modified Cholesky factor, E is the
  modification matrix which is diagonal and nonnegative, and P is the
  permutation matrix from pivoting.
  The input/output matrices are dense and stored in Matrix-Market format.
Arguments:
  INPUT_MATRIX_A is the file name of the input symmetric matrix A.
  OUTPUT_MATRIX_L is the file name of the output modified Cholesky factor L.
  The log will print out the permutation corresponding to P, as well as the
  the diagonal of P*E*P^T (the modification).
Optional parameters:
  -P=MATRIX_FILE      save the permutation matrix P to MATRIX_FILE
  -E=MATRIX_FILE      save the diagonal modification matrix E to MATRIX_FILE
  -help               display this help and exit
  -se90               for the SE90 algorithm
  -se99               for the SE99 algorithm
  -se1                for the SE-I algorithm
                      if "-se90",  "-se99", or "-se1" is specified,
                      then the following parameter specification is ignored
  -delta=VALUE        the modification tolerance delta
  -nopivot            for no pivoting in phase 2
                      (the default is Gershgorin circle pivoting in phase 2)
  -type2              for a type II algorithm (default is type I)
  -nondecreaseing     for the diagonal modifications being non-decreasing
                      (the default allows decreasing modifications)
  -relax=VALUE        relaxation factor, which must be positive and <= 1.0
                      (the default is not to use the relaxation strategy)
  -no_special_last    for not using the special treatment for the last
                      1-by-1 or 2-by-2 steps used in the SE99 algorithm
                      (the default invokes this special treatment)
