//----------------------------------------------------------------------------
// A mex wrapper for the modified LDL^T algorithm SE-I.
//
// Reference:
// Modified Cholesky Algorithms: A Catalog with New Approaches,
// Haw-ren Fang and Dianne O'Leary,
// Mathematical Programming, Series A, Vol. 115, No. 2, pp. 319--349, 2008
//----------------------------------------------------------------------------

#include "mex.h"
#include "../source/se.h"

extern
void mldl_mex_wrapper(const char *command,                                          // OCTAVE/MATLAB command
                      int (*mldl)(const int, double*, int *&, double *),            // routine to call
                      int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);  // mex I/O

// OCTAVE/MATLAB:  L        = se1(A);
//                [L, P]    = se1(A);
//                [L, P, E] = se1(A);
// ======
// input:
// A is a symmetric real matrix
// ======
// output:
// se1(A) factors input symmetric real matrix A as P*(A+E)*P'=L*L' by a modified LDL^T algorithm,
// where L is lower triangular
//       P is the permutation matrix from pivoting
//       E is the nonnegative diagonal modification matrix
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    const char command[] = "se1";
    mldl_mex_wrapper(command, &mchol_se1, nlhs, plhs, nrhs, prhs);
}
