//----------------------------------------------------------------------------
// A mex wrapper for modified LDL^T algorithms (i.e. the GMW and SE families).
//
// Reference:
// Modified Cholesky Algorithms: A Catalog with New Approaches,
// Haw-ren Fang and Dianne O'Leary,
// Mathematical Programming, Series A, Vol. 115, No. 2, pp. 319--349, 2008
//----------------------------------------------------------------------------

#include <stdio.h>
#include <string.h>
#include <limits>  // for std::numeric_limits<class T>::epsilon()
#include "mex.h"


void print_usage(const char *command);

// in all the comments, we assume command[] is "gmw81", but it can also be "gmw1", "gmw2", "se90", "se99", or "se1"

// OCTAVE/MATLAB:  L        = gmw81(A);
//                [L, P]    = gmw81(A);
//                [L, P, E] = gmw81(A);
// ======
// input:
// A is a symmetric real matrix
// ======
// output:
// gmw81(A) factors input symmetric real matrix A as P*(A+E)*P'=L*L' by a modified LDL^T algorithm,
// where L is lower triangular
//       P is the permutation matrix from pivoting
//       E is the nonnegative diagonal modification matrix
void mldl_mex_wrapper(const char *command,                                         // OCTAVE/MATLAB command
                      int (*mldl)(const int, double*, int *&, double *),           // routine to call
                      int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])  // mex I/O
{
    if (nrhs < 1) {
        print_usage(command);
        return;
    }

    // quick check of input matrix A
    if (mxGetNumberOfDimensions(prhs[0]) != 2) {
        mexPrintf("%s: input must be a matrix (i.e. 2-dimensional)!", command);
        mexErrMsgTxt("");
    }
    if (mxIsSparse(prhs[0])) {
        mexPrintf("%s: input cannot be a sparse matrix!", command);
        mexErrMsgTxt("");
    }
    if (mxIsComplex(prhs[0])) {
        mexPrintf("%s: input cannot be a complex matrix!", command);
        mexErrMsgTxt("");
    }
    if (!mxIsDouble(prhs[0])) {
        mexPrintf("%s: input must be double (i.e. cannot be integer, etc.)!", command);
        mexErrMsgTxt("");
    }

    // get the size of input matrix A
    size_t nrows = mxGetM(prhs[0]);
    size_t ncols = mxGetN(prhs[0]);
    // equivaliently in OCTAVE/MATLAB: [nrows, ncols] = size(A);
    if (nrows != ncols) {
        mexPrintf("%s: input matrix must be square!", command);
        mexErrMsgTxt("");
    }
    if (nrows<=0 || ncols<=0) {
        mexPrintf("%s: input matrix cannot be null!", command);
        mexErrMsgTxt("");
    }

    // get the address of input matrix A
    const double *A = mxGetPr(prhs[0]);

    // check whether A is symmetric
    for (int j=0; j<ncols; j++) {
        for (int i=j+1; i<nrows; i++) {
            if (A[j*nrows+i] != A[i*nrows+j]) {
                mexPrintf("%s: input matrix is not symmetric!", command);
                mexErrMsgTxt("");
            }
        }
    }

    // convert A to be in the compact format (only lower triangular part is stored)
    double *A2 = new double[nrows*(nrows+1)/2];
    double *a2 = A2;
    for (int j=0; j<ncols; j++) {
        for (int i=j; i<nrows; i++)
            *a2++ = A[j*nrows+i];
    }

    // perform the modified LDL^T algorithm
    int *perm = NULL;  // required memory perm[] will be allocated by new [] in mldl()
    double *modified = (nlhs==3) ? (new double[nrows]) : NULL;  // if nlhs==3, we also need to output the diagonal modification matrix E
    mldl(nrows, A2,
         perm, modified);
    // the algorithm is in-place; the modified Cholesky factor is stored in A2[]

    // set the output L,P,E from P*(A+E)*P' = L*L'
    // the L part (modified Cholesky factor)
    plhs[0] = mxCreateDoubleMatrix(nrows, ncols, mxREAL);  // nrows == ncols
    double *L = mxGetPr(plhs[0]);
    a2 = A2;
    for (int j=0; j<ncols; j++) {
        for (int i=j; i<nrows; i++) {
            L[i*nrows+j] = 0.0;
            L[j*nrows+i] = *a2++;
        }
    }
    // the P part (permutation matrix), stored as a sparse matrix (CSC format used by OCTAVE/MATLAB)
    if (nlhs >= 2) {
        plhs[1] = mxCreateSparse(nrows, nrows, nrows, mxREAL);
        double  *sp = mxGetPr(plhs[1]);
        mwIndex *ii = mxGetIr(plhs[1]);
        mwIndex *jj = mxGetJc(plhs[1]);
        for (int i=0; i<nrows; i++) {
            sp[i] = 1.0;      // value
            ii[perm[i]] = i;  // row index
            jj[i] = i;        // column pointer
        }
        jj[nrows] = nrows;    // terminator
    }
    // the E part (diagonal modification matrix)
    if (nlhs >= 3) {
        plhs[2] = mxCreateSparse(nrows, nrows, nrows, mxREAL);
        double  *se = mxGetPr(plhs[2]);
        mwIndex *ii = mxGetIr(plhs[2]);
        mwIndex *jj = mxGetJc(plhs[2]);
        for (int i=0; i<nrows; i++) {
            se[i] = modified[i];  // value
            ii[i] = i;            // row index
            jj[i] = i;            // column pointer
        }
        jj[nrows] = nrows;        // terminator
    }

    // free memory
    delete [] A2;
    delete [] perm;
    delete [] modified;  // it is safe to delete NULL
}


void print_usage(const char *command)
{
    mexPrintf("Usage:  L        = %s(A);\n", command);
    mexPrintf("       [L, P]    = %s(A);\n", command);
    mexPrintf("       [L, P, E] = %s(A);\n", command);
    mexPrintf("Input: A is a symmetric real matrix.\n");
    mexPrintf("Output: The factorization by the %s algorithm P*(A+E)*P'=L*L', where\n", command);
    mexPrintf("        L is lower triangular,\n");
    mexPrintf("        P is the permutation matrix from pivoting, and\n");
    mexPrintf("        E is the nonnegative diagonal modification matrix.\n");
}
