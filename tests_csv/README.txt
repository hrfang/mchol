This set of codes (gmw81, gmw1, gmw2, se90, se99, and se1) compute Cholesky
factorizations of real symmetric matrices, modified if necessary to make
them positive definite. They differ in how they compute the modification.

Given an input matrix A, each produces output matrices L, P, and E so that 
            P * (A + E) * P' = L * L'
where P is a permutation matrix, L is lower triangular, and E is the
modification.


======

The tests here require the 6 binary executable mex files for Octave or
Matlab,
------
../mex/gmw81.mex
../mex/gmw1.mex
../mex/gmw2.mex
../mex/se90.mex
../mex/se99.mex
../mex/se1.mex
------
which should be available after successful installation by executing "make"
under Octave or Matlab. For Matlab 64bit machines, the mex files may be
have extension "mexa64" instead of "mex" (i.e. "gmw81.mexa64" instead of
"gmw81.mex", etc.). The compiled binary executable files in "../drivers/"
are not required.

In contrast, the tests in "../tests_mtx/", which parallel the tests here,
require the binary executable files in "../drivers/" but not the mex files.

The codes in the tests here are verified working in both Matlab R2015a (trial
version) and Octave 3.8.1.

There are 1+33 test matrices included in the tests. These matrices are used
for testing the modified Cholesky algorithms in the SE99 paper as well as
our paper. The SE99 paper refers to:
------
A Revised Modified Cholesky Factorization Algorithm,
Robert B. Schnabel and Elizabeth Eskow,   
SIAM J. Optim., Vol. 9 No. 4, 1135-1148, 1999.


======

The two m-files to read and write matrices in matrix-market format are from
here:
http://math.nist.gov/MatrixMarket/mmio/matlab/mmread.m
http://math.nist.gov/MatrixMarket/mmio/matlab/mmwrite.m

Note that the logical and "&" and or "|" are replaced by "&&" and "||" to
get rid of the warnings in Octave. (In Octave, "&", "|" are bitwise "and",
"or", respectively.)

In the following tests, we need "mmread.m" to read matrix files in
matrix-market format. We don't need "mmwrite.m", which is included here for
convenience.


======

"csv_list.txt" lists all matrices for testing.

"test_gmw_fac.m" applies GMW81, GMW-I, GMW-II algorithms to the test
matrices listed in "csv_list.txt" (1+33 matrices in total).

For example, "../data/csv/first.csv" is in "csv_list.txt".
"test_gmw_fac.m" executes (via the mex file):
------
A = csvread(fileA);  % fileA is "../data/csv/first.csv"
[L, P, E] = gmw81(A);
------
The factorization is in the form: P*(A+E)*P' = L*L', where
L is the modified Cholesky factor stored in "FIRST_gmw81_L.csv",
P is the permutation matrix from pivoting stored in "FIRST_gmw81_P.csv", and
E is the diagonal modification matrix stored in "FIRST_gmw81_E.csv".

In addition to the GMW81 algorithm, GMW-I and GMW-II are also applied with
the substring "gmw81" replaced by "gmw1" and "gmw2", respectively.

The three algorithms repeated for the other 33 matrices.

The log file "gmw_result/log/test_gmw_fac.log" is from (in Octave) executing
------
diary test_gmw_fac.log
test_gmw_fac
diary off
------


======

"test_gmw_mod.m" checks the modification results from "test_gmw_fac.m". For
each combination of algorithm (GMW81, GMW-I, GMW-II) and test matrix, it
reports:
1) r2 = norm(E)/min{|lambda(A)|}, where norm(E) is the 2-norm of E, and
   min{|lambda(A)|} is least modification measured in 2-norm to make A+E
   positive semi-definite. An exception is that when A is already positive
   semi-definite, r2 is set as norm(E).
2) zeta = floor(kappa_2(A+E)), where kappa_2(A+E) is the condition number of
   A+E (2-norm based).

With the output matrices from executing "test_gmw_fac" in Octave, the log
file "gmw_result/log/test_gmw_mod.log" is from (in Octave) executing
------
diary test_gmw_mod.log
test_gmw_mod
diary off
------

To browse the result of one particular method (e.g. GMW81), one can execute
the Linux/Unix command:
------
grep "proc\|gmw81" test_gmw_mod.log > gmw81_mod.log
------

The procedure can be repeated for the results of SE90, SE99, SE-I by the 2
m-files, "test_se_fac.m" and "test_se_mod.m".

The sample results by Octave are included in "./result_octave/", and the
sample results by Matlab are included in "./result_matlab/".


======

The results match those reported in our paper (Tables 5,6,7):
Modified Cholesky Algorithms: A Catalog with New Approaches,
Haw-ren Fang and Dianne O'Leary,
Mathematical Programming, Series A, Vol. 115, No. 2, pp. 319--349, 2008
