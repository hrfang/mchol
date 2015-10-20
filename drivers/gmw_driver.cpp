//----------------------------------------------------------------------------
// A driver for modified Cholesky algorithms in the GMW family.
//
// Reference:
// Modified Cholesky Algorithms: A Catalog with New Approaches,
// Haw-ren Fang and Dianne O'Leary,
// Mathematical Programming, Series A, Vol. 115, No. 2, pp. 319--349, 2008
//----------------------------------------------------------------------------

#include <stdlib.h>  // for exit, atof, atoi, atol, srand, rand, etc.
#include <string.h>  // for memcpy, strcmp, strncmp, strlen, etc.
#include <time.h>    // for time_t, time, clock_t, clock, CLOCKS_PER_SEC, etc.
#include <iostream>  // for cout, cerr, endl, etc. (under namespace std)
#include <iomanip>   // for setw, setprecision, etc. (under namespace std)

#include "mtx_io.h"
#include "gmw.h"

using std::cout;
using std::cerr;
using std::endl;

using std::setw;
using std::setprecision;
using std::scientific;

void print_usage_and_exit(const char *cmd);

#define REAL double

int main(int argc, char *argv[]) {
    char *cmd = argv[0];
    if (argc <= 2 ||                             // argc <= 1  means no input matrix is specified
        !strcmp(argv[1], "-?") || !strcmp(argv[1], "-h") || !strcmp(argv[1], "-help") || !strcmp(argv[1], "--help"))
        print_usage_and_exit(cmd);               // print usage and exit

    // declare variables for input parameters
    REAL delta = 0.0;
    int pivot_method = 0;     // 0  for no pivoting;
                              // 1  for pivoting by maximum diagonal value;
                              // 2  for pivoting by maximum diagonal magnitude
    bool is_type1 = true;
    bool nondecreasing = false;
    bool is_2phase = false;
    REAL relax_factor = 0.0;  // default 0.0 means not to use the relaxation strategy
                              // relaxation strategy is the key factor that the SE99 algorithm improves the stability of the SE90 algorithm
                              // so it is generally good to set the relaxation factor, which must be >0.0 and <=1.0
                              // the SE99 algorithm uses relax_factor=1.0
                              // for GMW, relax_factor=0.75 achieves better empirical result for the 33 matrices
    bool special_last = false;
    int  algo = -1;           // -1 for unspecified
                              // 0  for GMW81  ("-gmw81")
                              // 1  for GMW-I  ("-gmw1")
                              // 2  for GMW-II ("-gmw2")

    // parse input parameters
    const char *input_A_file  = NULL;
    const char *output_L_file = NULL;
    const char *output_P_file = NULL;
    const char *output_E_file = NULL;
    for (int i=1; i<argc; i++) {
        // particular members of the GMW family: GMW81, GMW-I, GMW-II
        if (!strcmp(argv[i], "-gmw81")) {
            if (algo != -1) {
                cout << "Error: algo (\"-algo=\") is specified more than once!" << endl;
                exit(1);
            }
            algo = 0;
        }
        else if (!strcmp(argv[i], "-gmw1")) {
            if (algo != -1) {
                cout << "Error: algo (\"-algo=\") is specified more than once!" << endl;
                exit(1);
            }
            algo = 1;
        }
        else if (!strcmp(argv[i], "-gmw2")) {
            if (algo != -1) {
                cout << "Error: algo (\"-algo=\") is specified more than once!" << endl;
                exit(1);
            }
            algo = 2;
        }
        // modification delta
        else if (!strncmp(argv[i], "-delta=", strlen("-delta="))) {
            if (delta != 0) {
                cout << "Error: delta (\"-delta=\") is specified more than once!" << endl;
                exit(1);
            }
            delta = atof(argv[i] + strlen("-delta="));  // delta
            if (delta <= 0.0) {
                cout << "Error: delta (\"-delta=\") must be positive!" << endl;
                exit(1);
            }
        }
        // pivoting method
        else if (!strncmp(argv[i], "-pivot=", strlen("-pivot="))) {
            if (pivot_method != 0) {
                cout << "Error: algo (\"-pivot=\") is specified more than once!" << endl;
                exit(1);
            }
            pivot_method = atoi(argv[i] + strlen("-pivot="));  // pivoting method
            if (pivot_method < 1 || pivot_method > 2) {
                cout << "Error: the pivoting method index (\"-pivot=\") must be 1 or 2!" << endl;
                exit(1);
            }
        }
        // type-II
        else if (!strcmp(argv[i], "-type2")) {
            if (!is_type1) {
                cout << "Error: type-II algorithm (\"-type2\") is specified more than once!" << endl;
                exit(1);
            }
            is_type1 = false;
        }
        else if (!strcmp(argv[i], "-nondecreasing")) {
            if (nondecreasing) {
                cout << "Error: nondecreasing strategy (\"-nondecreasning\") is specified more than once!" << endl;
                exit(1);
            }
            nondecreasing = true;
        }
        // relaxation factor
        else if (!strncmp(argv[i], "-relax=", strlen("-relax="))) {
            if (relax_factor != 0.0) {
                cout << "Error: the relaxation factor (\"-relax=\") is specified more than once!" << endl;
                exit(1);
            }
            relax_factor = atof(argv[i] + strlen("-relax="));  // relaxation factor of the 2-phase strategy
            if (relax_factor <= 0.0 || relax_factor > 1.0) {
                cout << "Error: the relaxation factor (\"-relax=\") must be positive and <= 1!" << endl;
                exit(1);
            }
        }
        // special treatment for the last 1-by-1 or 2-by-2 Schur complement
        else if (!strcmp(argv[i], "-special_last")) {
            if (special_last) {
                cout << "Error: the special treatment for the last 1-by-1 or 2-by-2 Schur complement (\"-special_last\") is specified more than once!" << endl;
                exit(1);
            }
            special_last = true;
        }
        // permutation matrix P in the factorization P*(A+E)*P^T = L*L^T
        else if (!strncmp(argv[i], "-P=", strlen("-P="))) {
            if (output_P_file) {
                cout << "Error: the output permutation matrix P (\"-P=\") is specified more than once!" << endl;
                exit(1);
            }
            output_P_file = argv[i] + strlen("-P=");
        }
        // diagonal modification matrix E in the factorization P*(A+E)*P^T = L*L^T
        else if (!strncmp(argv[i], "-E=", strlen("-E="))) {
            if (output_E_file) {
                cout << "Error: the output diagonal matrix E (\"-E=\") is specified more than once!" << endl;
                exit(1);
            }
            output_E_file = argv[i] + strlen("-E=");
        }
        // invalid argument
        else if (!strncmp(argv[i], "-", strlen("-"))) {
            cerr << "Error: the argument \"" << argv[i] << "\" is not recognized!" << endl;
            exit(1);
        }
        else if (!input_A_file) {
            input_A_file = argv[i];
        }
        else if (!output_L_file) {
            output_L_file = argv[i];
        }
        else {
            cerr << "Error: more than 2 arguments as mtx file names are specified (exactly 2 mtx file names are needed, one for input A and the other for output L)!" << endl;
            exit(1);
        }
    }
    if (!input_A_file) {
        cerr << "Error: input matrix A file is not specified!" << endl;
        exit(1);
    }
    if (!output_L_file) {
        cerr << "Error: output matrix L file is not specified!" << endl;
        exit(1);
    }

    // open the matrix file and read the matrix
    cout << "------" << endl;
    cout << "mtx file of input  A: " << input_A_file  << endl;
    cout << "mtx file of output L: " << output_L_file << endl;
    if (output_P_file)
        cout << "mtx file of output P: " << output_P_file << endl;
    int nrows, ncols;
    char symm;
    REAL *sa = NULL;
    const bool convert = true;
    int info = matrix_market_matrix_read(input_A_file,            // input
                                         nrows, ncols, sa, symm,  // output
                                         convert);                // optional parameter
    if (info) {
        cerr << "Error: failed to read input mtx file \"" << input_A_file << "\" (info=" << info << ")!" << endl;
        exit(1);
    }
    cout << "nrows=" << nrows << ", ncols=" << ncols << ", symm='" << symm << "'" << endl;

    // make sure the matrix is symmetric
    if (symm!='s' && symm!='S') {
        cerr << "Error: input matrix is not symmetric!" << endl;
        exit(1);
    }

    int *perm = new int[nrows];
    for (int i=0; i<nrows; i++)
        perm[i] = i;
    REAL *modified = new REAL[nrows];  // for storing modified diagonal elements ("diagonal" is after pivoting)

    if (algo == 0) {
        // GMW81 algorithm
        info = mchol_gmw81(nrows, sa,
                          perm, modified);
    }
    else if (algo == 1) {
        // GMW-I algorithm
        info = mchol_gmw1(nrows, sa,
                          perm, modified);
    }
    else if (algo == 2) {
        // GMW-II algorithm
        info = mchol_gmw2(nrows, sa,
                          perm, modified);
    }
    else {
        // algo==-1, a generic member in the GMW family
        info = mchol_gmw(nrows, sa,
                         perm, modified,
                         delta,
                         pivot_method,
                         is_type1, nondecreasing,
                         is_2phase, relax_factor,
                         special_last);
    }

    if (info) {
        cerr << "Error: modified Cholesky factorization failed (info=" << info << ")!" << endl;
        exit(1);
    }

    // print modified[]
    cerr << "modified:";
    for (int i=0; i<nrows; i++)
        cerr << " " << modified[i];
    cerr << endl;

    // print perm[]
    cerr << "perm:";
    for (int i=0; i<nrows; i++)
        cerr << " " << perm[i]+1;
    cerr << endl;

    info = matrix_market_matrix_write(nrows, ncols, sa, 'l',  // 'l' means that the matrix stored in sa[] is lower triangular
                                      output_L_file);
    if (info) {
        cerr << "Error: failed to write output L to mtx file \"" << output_L_file << "\" (info=" << info << ")!" << endl;
        exit(1);
    }
    if (output_P_file) {
        // allocate memory for sp[] and zero initialization
        int sz = nrows*ncols;
        REAL *sp = new REAL[sz];
        REAL *sp0 = sp;
        while (sz--)
            *sp0++ = 0.0;
        // now sp[i]==0.0 for i=0,...,nrows*ncols-1
        for (int i=0; i<nrows; i++)
            sp[perm[i]*nrows+i] = 1.0;
        info = matrix_market_matrix_write(nrows, ncols, sp, 'g',  // 'g' for that P is saved as a general matrix
                                          output_P_file);
        if (info) {
            cerr << "Error: failed to write output permutation matrix P to mtx file \"" << output_P_file << "\" (info=" << info << ")!" << endl;
            exit(1);
        }
        delete [] sp;
    }
    if (output_E_file) {
        info = matrix_market_diagonal_matrix_write(nrows, ncols, modified, output_E_file);
        if (info) {
            cerr << "Error: failed to write diagonal modification matrix E to mtx file \"" << output_E_file << "\" (info=" << info << ")!" << endl;
            exit(1);
        }
    }

    // free memory
    delete [] sa;
    delete [] modified;
    delete [] perm;

    return 0;
}

void print_usage_and_exit(const char *cmd) {
    cout << "Usage: " << cmd << " INPUT_MATRIX_A OUTPUT_MATRIX_L [OPTION] [OPTION] ..." << endl;
    cout << "  Modified Cholesky algorithms in the GMW family." << endl;
    cout << "  Given the input symmtric matrix A, the factorization is in the form" << endl;
    cout << "                P*(A+E)*P^T = L*L^T," << endl;
    cout << "  where L is the lower triangular modified Cholesky factor, E is the" << endl;
    cout << "  modification matrix which is diagonal and nonnegative, and P is the" << endl;
    cout << "  permutation matrix from pivoting." << endl;
    cout << "  The input/output matrices are dense and stored in Matrix-Market format." << endl;
    cout << "Arguments:" << endl;
    cout << "  INPUT_MATRIX_A is the file name of the input symmetric matrix A." << endl;
    cout << "  OUTPUT_MATRIX_L is the file name of the output modified Cholesky factor L." << endl;
    cout << "  The log will print out the permutation corresponding to P, as well as the" << endl;
    cout << "  the diagonal of P*E*P^T (the modification)." << endl;
    cout << "Optional parameters:" << endl;
    cout << "  -P=MATRIX_FILE      save the permutation matrix P to file MATRIX_FILE" << endl;
    cout << "  -E=MATRIX_FILE      save the diagonal modification matrix E to MATRIX_FILE" << endl;
    cout << "  -help               display this help and exit" << endl;
    cout << "  -gmw81              for the GMW81  algorithm" << endl;
    cout << "  -gmw1               for the GMW-I  algorithm" << endl;
    cout << "  -gmw2               for the GMW-II algorithm" << endl;
    cout << "                      if \"-gmw81\",  \"-gmw1\", or \"-gmw2\" is specified," << endl;
    cout << "                      then the following parameter specification is ignored" << endl;
    cout << "  -delta=VALUE        the modification tolerance delta" << endl;
    cout << "  -pivot=INDEX        the pivoting method (default is no pivoting)" << endl;
    cout << "                      1 for pivoting by maximum diagonal element," << endl;
    cout << "                      2 for pivoting by maximum diagonal magnitude" << endl;
    cout << "  -type2              for a type II algorithm (default is type I)" << endl;
    cout << "  -nondecreaseing     for the diagonal modifications being non-decreasing" << endl;
    cout << "                      (the default allows decreasing modifications)" << endl;
    cout << "  -relax=VALUE        relaxation factor, which must be positive and <= 1.0" << endl;
    cout << "                      (the default is not to use the relaxation strategy)" << endl;
    cout << "  -special_last       for invoking the special treatment for the last 1-by-1" << endl;
    cout << "                      or 2-by-2 steps used in the SE99 algorithm" << endl;
    cout << "                      (the default is not to use this special treatment)" << endl;
    cout << endl;

    exit(1);
}
