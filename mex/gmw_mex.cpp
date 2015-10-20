//----------------------------------------------------------------------------
// A mex wrapper for the modified Cholesky algorithms in the GMW family.
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
#include "../source/gmw.h"

void print_usage();
void parse_gmw_options(const mxArray *mx_opts,  // input
                       double &delta,           // the rest are output parameters
                       int &pivot_method,
                       bool &is_type1,
                       bool &nondecreasing,
                       bool &is_2phase,
                       double &relax_factor,
                       bool &special_last);

// OCTAVE/MATLAB:  L        = gmw(A, opts);
//                [L, P]    = gmw(A, opts);
//                [L, P, E] = gmw(A, opts);
// ======
// input:
// A is a symmetric real matrix
// ======
// output:
// gmw(A, opts) factors input symmetric real matrix A as P*(A+E)*P'=L*L' by the GMW algorithm,
// where L is lower triangular
//       P is the permutation matrix from pivoting
//       E is the nonnegative diagonal modification matrix
// ======
// options (algorithm arguments):
// opts (struct) is a collection of algorithm parameters; see the header description of parse_gmw_options() for detail
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs < 2) {
        print_usage();
        return;
    }

    // quick check of input matrix A
    if (mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgTxt("gmw: input must be a matrix (i.e. 2-dimensional)!");
    if (mxIsSparse(prhs[0]))
        mexErrMsgTxt("gmw: input cannot be a sparse matrix!");
    if (mxIsComplex(prhs[0]))
        mexErrMsgTxt("gmw: input cannot be a complex matrix!");
    if (!mxIsDouble(prhs[0]))
        mexErrMsgTxt("gmw: input must be double (i.e. cannot be integer, etc.)!");

    // get the size of input matrix A
    size_t nrows = mxGetM(prhs[0]);
    size_t ncols = mxGetN(prhs[0]);
    // equivaliently in OCTAVE/MATLAB: [nrows, ncols] = size(A);
    if (nrows != ncols)
        mexErrMsgTxt("gmw: input matrix must be square!");
    if (nrows<=0 || ncols<=0)
        mexErrMsgTxt("gmw: input matrix cannot be null!");

    // get the address of input matrix A
    const double *A = mxGetPr(prhs[0]);

    // check whether A is symmetric
    for (int j=0; j<ncols; j++) {
        for (int i=j+1; i<nrows; i++) {
            if (A[j*nrows+i] != A[i*nrows+j])
                mexErrMsgTxt("gmw: input matrix is not symmetric!");
        }
    }

    // convert A to be in the compact format (only lower triangular part is stored)
    double *A2 = new double[nrows*(nrows+1)/2];
    double *a2 = A2;
    for (int j=0; j<ncols; j++) {
        for (int i=j; i<nrows; i++)
            *a2++ = A[j*nrows+i];
    }

    // get the GMW algorithm parameters from the 2nd input argument
    double delta;
    int pivot_method;
    bool is_type1;
    bool nondecreasing;
    bool is_2phase;
    double relax_factor;
    bool special_last;
    parse_gmw_options(prhs[1],
                      delta, pivot_method, is_type1, nondecreasing, is_2phase, relax_factor, special_last);

    // perform the GMW algorithm
    int *perm = NULL;  // required memory perm[] will be allocated by new [] in mchol_gmw()
    double *modified = (nlhs==3) ? (new double[nrows]) : NULL;  // if nlhs==3, we also need to output the diagonal modification matrix E
    mchol_gmw(nrows, A2,
              perm, modified,
              delta, pivot_method, is_type1, nondecreasing, is_2phase, relax_factor, special_last);
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


void print_usage()
{
    mexPrintf("Usage:  L        = gmw(A, opts);\n");
    mexPrintf("       [L, P]    = gmw(A, opts);\n");
    mexPrintf("       [L, P, E] = gmw(A, opts);\n");
    mexPrintf("Input: A is a symmetric real matrix.\n");
    mexPrintf("Output: The factorization by the GMW algorithm P*(A+E)*P'=L*L', where\n");
    mexPrintf("        L is lower triangular,\n");
    mexPrintf("        P is the permutation matrix from pivoting, and\n");
    mexPrintf("        E is the nonnegative diagonal modification matrix.\n");
    mexPrintf("Options (algorithm arguments):\n");
    mexPrintf("    opts.delta: the modification tolerance, which must be positive\n");
    mexPrintf("                the default delta is eps (machine epsilon)\n");
    mexPrintf("    opts.pivot_method: 0 for no pivoting (default)\n");
    mexPrintf("                       1 for pivoting by maximum diagonal element\n");
    mexPrintf("                       2 for pivoting by maximum diagonal magnitude\n");
    mexPrintf("    opts.is_type1: true  for type 1 modification algorithm (default)\n");
    mexPrintf("                   false for type 2 modification algorithm\n");
    mexPrintf("    opts.nondecreaseing: true  for the nondecreasing modification\n");
    mexPrintf("                         false for not enforcing it (default)\n");
    mexPrintf("    opts.is_2phase: true  for the 2-phase strategy\n");
    mexPrintf("                    false for not applying it (default)\n");
    mexPrintf("    opts.relax_factor: the relaxation factor in phase 1\n");
    mexPrintf("                       the value must be in (0,1]\n");
    mexPrintf("                       effective only if opts.is_2phase==true, in which case\n");
    mexPrintf("       if relax_factor is not specified, the 2-phase strategy is not relaxed\n");
    mexPrintf("    opts.special_last: true  for the SE special treatment in the last step\n");
    mexPrintf("                       false for not employing it (default)\n");
}


// convert an input argument mxArray "in" of OCTAVE/MATLAB to "val" as a double precision number
// the converted value "val" is double, even if the input "in" is float or int
// on success, return 0; otherwise, the input "in" is not numeric and return -1
static
int mxArray_to_double(const mxArray *in,  // input
                      double &val,        // output
                      int offset=0)       // optional offset
{
    // return -1 if input in is not numeric
    if (!mxIsNumeric(in))
        return -1;

    void *dat = mxGetData(in);
    switch (mxGetClassID(in)) {
        case mxDOUBLE_CLASS:
            val = static_cast<double>(*((double *)(dat)+offset));
            break;
        case mxSINGLE_CLASS:
            val = static_cast<double>(*((float *)(dat)+offset));
            break;
        case mxINT8_CLASS:
            val = static_cast<double>(*((int8_t *)(dat)+offset));
            break;
        case mxUINT8_CLASS:
            // val = static_cast<double>(*((uint8_t *)(dat)+offset));  // "uint8_t" not recognized by OCTAVE compiler mkoctfile ver 3.2.4
            val = static_cast<double>(*((unsigned char *)(dat)+offset));
            break;
        case mxINT16_CLASS:
            val = static_cast<double>(*((int16_t *)(dat)+offset));
            break;
        case mxUINT16_CLASS:
            // val = static_cast<double>(*((uint16_t *)(dat)+offset));  // "uint16_t" not recognized by OCTAVE compiler mkoctfile ver 3.2.4
            val = static_cast<double>(*((unsigned short *)(dat)+offset));
            break;
        case mxINT32_CLASS:
            val = static_cast<double>(*((int32_t *)(dat)+offset));
            break;
        case mxUINT32_CLASS:
            // return static_cast<double>(*((uint32_t *)(dat)+offset));  // "uint32_t" not recognized by OCTAVE compiler mkoctfile ver 3.2.4
            val = static_cast<double>(*((unsigned *)(dat)+offset));
            break;
        default:
            // the remains are not numeric (cell, string, or structure) and has been handled in "if (!mxIsNumeric(ma)" above
            // should not end in here unless new types of mxArray are introduced
            return -1;
    }

    return 0;  // success
}

// mx_opts is the input (struct) mxArray
// the rest arguments are the output algorithm parameters from mx_opts
//    delta is the modification tolerance, which must be positive
//       if option delta is not present, it will be set as eps
//    pivot_method = 0, for no pivoting (default)
//                 = 1, for pivoting by maximum diagonal element
//                 = 2, for pivoting by maximum diagonal magnitude
//    is_type1 indicates whether the modification is of type 1 (default) or type 2
//    nondecreaseing indicates whether the nondecreasing modification strategy is enforced or not (default is not)
//    is_2phase indicates whether the 2-phase strategy is applied or not (default is not)
//       if 2-phase strategy is applied, then phase 1 is always pivoted by the maximum diagonal element
//    relax_factor is effective only if the 2-phase strategy is applied (i.e. is_2phase==true)
//       relax_factor must be in (0,1] as the relaxation factor (i.e. mu in our 2008 modified Cholesky paper)
//       if relax_factor is not specified, it means that the 2-phase strategy is not relaxed
//    special_last indicates whether the SE special treatment for the last 1-by-1 or 2-by-2 Schur complement is applied
// if mx_opts contains an invalid parameter, an error message will be printed via mexErrMsgTxt()
void parse_gmw_options(const mxArray *mx_opts,
                       double &delta,
                       int &pivot_method,
                       bool &is_type1,
                       bool &nondecreasing,
                       bool &is_2phase,
                       double &relax_factor,
                       bool &special_last)
{
    // get input arguments
    int num_fields = mxGetNumberOfFields(mx_opts);
    if (mxGetNumberOfElements(mx_opts) != 1) {
        mexErrMsgTxt("gmw: option error, each field should have one and only one element!");
    }

    // set default algorithm parameters
    delta = std::numeric_limits<double>::epsilon();  // eps (machine epsilon)
    pivot_method = 0;
    is_type1 = true;
    nondecreasing = false;
    is_2phase = false;
    relax_factor = 0.0;
    special_last = false;

    // parse the options
    double val;
    for (int ifield=0; ifield<num_fields; ifield++) {
        size_t nrows, ncols;
        mxArray *field_data = mxGetFieldByNumber(mx_opts, 0, ifield);
        // if the field is empty, ignore
        if ((nrows=mxGetM(field_data))==0 || (ncols=mxGetN(field_data))==0)
            continue;
        // if the field is not a single entry, print an error message
        if (nrows!=1 || ncols!=1)
            mexErrMsgTxt("gmw: all options should be single numbers (cannot be a matrix)!");
        // if the field is not numeric, print an error message
        const char *field_name = mxGetFieldNameByNumber(mx_opts, ifield);
        if (mxArray_to_double(field_data, val) == -1)
            mexErrMsgTxt("gmw: all options must be numeric (cannot be a cell or structure)!");

        if (!strcmp(field_name, "delta")) {
            if (val < 0.0)
                mexErrMsgTxt("gmw: the option delta must be positive!");
            delta = val;
        }
        else if (!strcmp(field_name, "pivot_method")) {
            if (val!=0.0 && val!=1.0 && val!=2.0)
                mexErrMsgTxt("gmw: the option pivot_method must be 0 (for no pivoting), 1 (for pivoting by max diagonal element), or 2 (for pivoting by max diagonal magnitude)!");
            pivot_method = static_cast<int>(val);
        }
        else if (!strcmp(field_name, "is_type1")) {
            is_type1 = (val!=0.0);
        }
        else if (!strcmp(field_name, "nondecreasing")) {
            nondecreasing = (val!=0.0);
        }
        else if (!strcmp(field_name, "is_2phase")) {
            is_2phase = (val!=0.0);
        }
        else if (!strcmp(field_name, "relax_factor")) {
            if (val<=0.0 || val>1.0)
                mexErrMsgTxt("gmw: the option relax_factor (relaxation factor) must be in (0,1]!");
            relax_factor = val;
        }
        else if (!strcmp(field_name, "special_last")) {
            special_last = (val!=0.0);
        }
        else {
            mexPrintf("gmw: invalid option '%s'!", field_name);
            mexErrMsgTxt("");
        }
    }
}
