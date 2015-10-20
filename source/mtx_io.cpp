//----------------------------------------------------------------------------
// This file contains routines to read or write dense matrices in matrix
// market format.
// This code comes with no guarantee or warranty of any kind.
//----------------------------------------------------------------------------

#include <stdio.h>   // for fopen, fclose, etc.
#include <stdlib.h>  // for exit, atof, atoi, atol, srand, rand, etc.
#include <string.h>  // for memcpy, strcmp, strncmp, strlen, etc.
#include <ctype.h>   // for tolower
#include <typeinfo>  // for typeid

#include "mtx_io.h"


// matrix_market_matrix_read() reads a real dense matrix from a matrix-market file file_name[],
// and store the result in the coordinate format in sa[], nrows, ncols,
// where nrows, ncols are the number of rows, the number of columns, respectively,
// and the entries are stored column-wisely in sa[]
// for a general real matrix A, the element A(i,j) is stored in sa[(j-1)*nrows+(i-1)]
// for a symmetric or skew-symmetric real matrix A, only the lower triangular part (including the diagonal) is stored
// symmetry information will also be recorded, with symm='g' for general, symm='s' for symmetric, and symm='k' for skew-symmetric
//    for the latter two cases ('s' and 'k'), each pair of off-diagonal elements A(i,j) and A(j,i) is recorded once in the output
// the return value:
//    0: a successful read
//    1: the input file is a sparse matrix
//    2: the input file is a complex (dense) matrix
//    3: the input file contains only sparsity pattern information
//   -1: file open error
//   -2: file is empty
//   -3: the first line (the header) of the file is empty
//   -4: invalid header
//   -5: invalid number of rows, number of columns, or for symmetric or skew-symmetric cases, nrows != ncols
//   -6: the matrix data is incomplete or invalid
// if sa==NULL is given, then the required memory will be allocated;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sa
// if convert==true, then a matrix stored as being general (symm=='g') will be checked whether it is symmetric
//    (symm=='s') or skew-symmetric (symm=='k'), and conversion is made whenever it is appropriate,
//    in which case symm=='S' or 'K' instead of 's' or 'k' is marked
template<class Real>
int matrix_market_matrix_read(const char file_name[],                         // input
                              int &nrows, int &ncols, Real *&sa, char &symm,  // output
                              const bool convert)                             // optional parameter
{
    // open the file to read
    FILE *mm_file = fopen(file_name, "r");
    if (mm_file == NULL) {
        // file open error!
        return -1;
    }

    char string[512];  // assume that there are at most 255 characters in each line of the header
    if (fgets(string, 512, mm_file) == NULL) {
        // file is empty!
        fclose(mm_file);
        return -2;
    }

    if (strlen(string) == 0) {
        // the first line in the input file is empty!
        fclose(mm_file);
        return -3;
    }

    if (*string != '%') {
        // assume that the file has flat coordinates and values
        symm = 'g';
    }
    else {
        char *token = strtok(string, " \n");  // `%%MatrixMarket' should be read
        if (token == NULL || strcmp(token, "%%MatrixMarket")) {
            // invalid header
            fclose(mm_file);
            return -4;
        }

        token = strtok(NULL, " \n");  // "matrix" should be read
        if (token == NULL) {
            // invalid header
            fclose(mm_file);
            return -4;
        }
        int len = static_cast<int>(strlen(token));
        for (int i=0; i<len; i++)
            token[i] = static_cast<char>(tolower(token[i]));
        if (strcmp(token, "matrix")) {
            // not "matrix"!?  invalid header
            fclose(mm_file);
            return -4;
        }

        token = strtok(NULL, " \n");  // "coordinate" should be read, for a sparse matrix
                                      // the other valid matrix-market tag is "array" for a general matrix
        if (token == NULL) {
            // invalid header
            fclose(mm_file);
            return -4;
        }
        len = static_cast<int>(strlen(token));
        for (int i=0; i<len; i++)
            token[i] = static_cast<char>(tolower(token[i]));
        if (!strcmp(token, "coordinate")) {
            // input matrix-market file stores a dense matrix; skip the rest
            fclose(mm_file);
            return 1;
        }
        if (strcmp(token, "array")) {
            // not "array"!?  invalid header
            fclose(mm_file);
            return -4;
        }

        token = strtok(NULL, " \n");  // "real" should be read, for a real-valued matrix, "integer" is also accepted
                                      // the other valid matrix-market tags are "complex" and "pattern"
        if (token == NULL) {
            // invalid header
            fclose(mm_file);
            return -4;
        }
        len = static_cast<int>(strlen(token));
        for (int i=0; i<len; i++)
            token[i] = static_cast<char>(tolower(token[i]));
        if (!strcmp(token, "complex")) {
            // the input matrix-market file stores a complex (sparse) matrix
            fclose(mm_file);
            return 2;
        }
        if (!strcmp(token, "pattern")) {
            // the input matrix-market file contains only the sparsity pattern information?
            fclose(mm_file);
            return 3;
        }
        if (strcmp(token, "real") && strcmp(token, "integer")) {
            // neither "real" nor "integer"!?  invalid header
            fclose(mm_file);
            return -4;
        }

        token = strtok(NULL, " \n");  // "general", "symmetric", "skew-symmetric" should be read, "hermitian" is also possible
                                      // "skew-hermitian" seems not a valid tag in the matrix-market format, but is allowed here
                                      // however, "hermitian" and "skew-hermitian" should not appear here, since we are reading a real-valued matrix
        if (token == NULL) {
            // invalid header
            fclose(mm_file);
            return -4;
        }
        len = static_cast<int>(strlen(token));
        for (int i=0; i<len; i++)
            token[i] = static_cast<char>(tolower(token[i]));
        if (!strcmp(token, "general")) {
            // input matrix-market file stores a general matrix
            symm = 'g';
        }
        else if (!strcmp(token, "symmetric") || !strcmp(token, "hermitian")) {
            // input matrix-market file stores a symmetric matrix
            symm = 's';
        }
        else if (!strcmp(token, "skew-symmetric") || !strcmp(token, "skew-hermitian")) {
            // input matrix-market file stores a skew-symmetric matrix
            symm = 'k';
        }
        else {
            // not a recognized tag;  invalid header
            fclose(mm_file);
            return -4;
        }

        // get rid of the rest comment lines
        do {
            if (fgets(string, 512, mm_file) == NULL) {
                // no information about nrows, ncols!?
                fclose(mm_file);
                return -5;
            }
        } while (strlen(string)==0 || *string=='\n' || *string=='%');
    }

    // now parse the line which contains nrows, ncols
    // read nrows
    char *token = strtok(string, " \n");
    if (token == NULL) {
        // no information about the number of rows!?
        fclose(mm_file);
        return -5;
    }
    nrows = atoi(token);
    if (nrows <= 0) {
        // number of rows is not positive!?
        fclose(mm_file);
        return -5;
    }

    // read ncols
    token = strtok(NULL, " \n");
    if (token == NULL) {
        // no information about the number of columns!?
        fclose(mm_file);
        return -5;
    }
    ncols = atoi(token);
    if (ncols <= 0) {
        // number of columns is not positive!?
        fclose(mm_file);
        return -5;
    }
    if ((symm=='s' || symm=='k') && (nrows!=ncols)) {
        // symmetric or skew-symmetric, but not square!?
        fclose(mm_file);
        return -5;
    }

    // allocate memory if necessary
    const int sz = (symm=='g') ? (nrows*ncols) : nrows*(nrows+1)/2;
    if (sa == NULL)
        sa = new Real[sz];

    // read the matrix
    Real *sa0 = sa;
    int sz0 = sz;
    while (sz0--) {
        int num_scanned = 0;
        if (typeid(Real) == typeid(double))
            num_scanned = fscanf(mm_file, "%lg", (double *)sa0);
        else  // assume "float"
            num_scanned = fscanf(mm_file, "%g",  (float *)sa0);
        sa0 ++;
        if (num_scanned != 1) {
            // incomplete data!?
            fclose(mm_file);
            return -6;
        }
    }
    if (symm == 'k') {
        // enforce zero diagonal of a skew-symmetric matrix
        sa0 = sa;
        for (int j=0; j<ncols; j++) {
            *sa0 = 0.0;
            sa0 += (nrows-j);
        }
    }

    // matrix conversion, if specified and appropriate
    if (convert && symm=='g' && nrows==ncols) {
        // check whether the matrix is symmetric
        sa0 = sa;
        symm = 'S';  // assume being symmetric; correct if it is not
        for (int j=0; j<ncols; j++) {
            const Real *sa1 = sa+j;  // address of A(1,j+1)
            for (int i=0; i<nrows; i++) {
                // *sa0 is A(i+1,j+1); *sa1 is A(j+1,i+1)
                if (*sa0 != *sa1) {
                    symm = 'g';
                    break;
                }
                sa0 ++;
                sa1 += nrows;
            }
        }
        if (symm == 'g') {
            // check whether the matrix is skew-symmetric
            sa0 = sa;
            symm = 'K';  // assume being symmetric; correct if it is not
            for (int j=0; j<ncols; j++) {
                const Real *sa1 = sa+j;  // address of A(1,j+1)
                for (int i=0; i<nrows; i++) {
                    // *sa0 is A(i+1,j+1); *sa1 is A(j+1,i+1)
                    if (*sa0 != -*sa1) {
                        symm = 'g';
                        break;
                    }
                    sa0 ++;
                    sa1 += nrows;
                }
            }
        }
        if (symm != 'g') {
            // symmetric or skew-symmtric, store only lower trianuglar part of the matrix
            sa0 = sa+nrows;  // for address of A(1,2)
            Real *sa1 = sa0;
            for (int j=1; j<ncols; j++) {
                sa1 += j;  // for address of A(j+1,j+1)
                int nr = nrows-j;
                while (nr--)
                    *sa0++ = *sa1 ++;
            }
        }
    }

    // a successful read, close the input file and return 0
    fclose(mm_file);
    return 0;
}

// matrix_market_matrix_write() writes a dense matrix A to a matrix-market file file_name[],
//    where A is stored column-wise in sa[], and
//    nrows, ncols are the number of rows, the number of columns, respectively, and
//    symm='g','s','k' for general, symmetric, skew-symmetric matrices, respectively
//    we also include the option symm='l', such that sa[] stores only the lower
//    triangular part of a lower triangular matrix columnwisely, but written to the
//    output file as a dense matrix (matrix-market format does not support the lower
//    triangular matrix)
// comment[] is a character string which will be printed in the output file file_name[] (if comment!=NULL)
// the return value:
//    0: a successful write
//    1: file open error
//    2: file close error
//    3: invalid symm (i.e. not 'g', 's', or 'k')
//    4: invalid nrows or ncols (i.e. not positive)
//    5: nrows!=ncols but symm='s' or 'k'
//    a negative integer: a number returned by fprintf() which signals an error
template<class Real>
int matrix_market_matrix_write(const int nrows, const int ncols, const Real *sa, const char symm,  // input
                               const char file_name[],                                             // output file name
                               const char comment[])                                               // optional comment
{
    // open the file to write
    FILE *mm_file = fopen(file_name, "w");
    if (mm_file == NULL) {
        // file open error!
        return 1;
    }

    // print the header
    int info;
    if (symm == 'g' || symm == 'G' || symm == 'l' || symm == 'L')
        info = fprintf(mm_file, "%%%%MatrixMarket matrix array real general\n");
    else if (symm == 's' || symm == 'S')
        info = fprintf(mm_file, "%%%%MatrixMarket matrix array real symmetric\n");
    else if (symm == 'k' || symm == 'K')
        info = fprintf(mm_file, "%%%%MatrixMarket matrix array real skew-symmetric\n");
    else {
        fclose(mm_file);
        return 3;
    }
    if (info < 0) {
        // file write error
        fclose(mm_file);
        return info;
    }

    // print the comment
    info = fprintf(mm_file, "%%-------------------------------------------------------------------------------\n");
    if (info >= 0) {
        if (comment == NULL)
            info = fprintf(mm_file, "%% matrix-market dense matrix file generated by matrix_market_matrix_write()\n");
        else
            info = fprintf(mm_file, "%% %s\n", comment);
        if (info >= 0)
            info = fprintf(mm_file, "%%-------------------------------------------------------------------------------\n");
    }
    if (info < 0) {
        // file write error
        fclose(mm_file);
        return info;
    }

    // print nrows (number of rows), ncols (number of columns)
    info = fprintf(mm_file, "%d %d\n", nrows, ncols);
    if (nrows<=0 || ncols<=0) {
        // invalid nrows or ncols
        fclose(mm_file);
        return 4;
    }
    if ((symm=='s' || symm=='S' || symm=='k' || symm=='K') && nrows!=ncols) {
        // a symmetric or skew-symmetric matrix should have nrows==ncols
        fclose(mm_file);
        return 5;
    }
    if (info < 0) {
        // file write error
        fclose(mm_file);
        return info;
    }

    // print the dense matrix to the output mtx file
    // matrix-market format is column-wise
    // so in text, if each line displays one column, it may look like transposed
    // (e.g. a lower triangular matrix may look like an upper triangular matrix)
    const Real *sa0 = sa;
    for (int j=0; j<ncols; j++) {
        const int nr = (symm=='g' || symm=='G') ? nrows : nrows-j;
        // if symm=='s','S','k','K','l','L' (symmetric, skew-symmetric, or lower-triangular),
        // then only the lower triangular part of the matrix is stored
        if (symm == 'l' || symm == 'L') {
            for (int i=0; i<j; i++) {
                info = fprintf(mm_file, " 0");
                if (info < 0) {
                    // file write error
                    fclose(mm_file);
                    return info;
                }
            }
        }
        for (int i=0; i<nr; i++) {
            if ((symm=='k' || symm=='K') && i==0)
                info = fprintf(mm_file, " 0");  // zero diagonal of a skew-symmetric matrix
            else if (typeid(Real) == typeid(double))
                info = fprintf(mm_file, " %.16lg", *sa0);
            else  // assume "float"
                info = fprintf(mm_file, " %g", *sa0);
            sa0 ++;
            if (info < 0) {
                // file write error
                fclose(mm_file);
                return info;
            }
        }
        info = fprintf(mm_file, "\n");
        if (info < 0) {
            // file write error
            fclose(mm_file);
            return info;
        }
    }

    // done writing; now close the file
    if (fclose(mm_file) == EOF) {
        // close file error
        return 2;
    }

    // a successful write, return 0
    return 0;
}


// instantiation by double

template
int matrix_market_matrix_read(const char file_name[],                           // input
                              int &nrows, int &ncols, double *&sa, char &symm,  // output
                              const bool convert);

template
int matrix_market_matrix_write(const int nrows, const int ncols, const double *sa, const char symm,  // input
                               const char file_name[],                                               // output file name
                               const char comment[]);                                                // optional comment

// instantiation by float

template
int matrix_market_matrix_read(const char file_name[],                          // input
                              int &nrows, int &ncols, float *&sa, char &symm,  // output
                              const bool convert);

template
int matrix_market_matrix_write(const int nrows, const int ncols, const float *sa, const char symm,  // input
                               const char file_name[],                                              // output file name
                               const char comment[]);                                               // optional comment


// matrix_market_diagonal_matrix_write() writes a diagonal matrix D to file_name[] in matrix-market sparse matrix format
//    where the diagonal of D is stored in sd[], and
//    nrows, ncols are the number of rows, the number of columns, respectively, and
// the return value:
//    0: a successful write
//    1: file open error
//    2: file close error
//    3: invalid nrc (i.e. not positive)
//    a negative integer: a number returned by fprintf() which signals an error
template<class Real>
int matrix_market_diagonal_matrix_write(const int nrows, const int ncols, const Real *sd,  // input
                                        const char file_name[],                            // output file name
                                        const char comment[])                              // optional comment
{
    // open the file to write
    FILE *mm_file = fopen(file_name, "w");
    if (mm_file == NULL) {
        // file open error!
        return 1;
    }

    // print the header
    int info = fprintf(mm_file, "%%%%MatrixMarket matrix coordinate real general\n");
    if (info < 0) {
        // file write error
        fclose(mm_file);
        return info;
    }

    // print the comment
    info = fprintf(mm_file, "%%-------------------------------------------------------------------------------\n");
    if (info >= 0) {
        if (comment == NULL)
            info = fprintf(mm_file, "%% matrix-market dense matrix file generated by matrix_market_diagonal_matrix_write()\n");
        else
            info = fprintf(mm_file, "%% %s\n", comment);
        if (info >= 0)
            info = fprintf(mm_file, "%%-------------------------------------------------------------------------------\n");
    }
    if (info < 0) {
        // file write error
        fclose(mm_file);
        return info;
    }

    const int nrc = (nrows<=ncols) ? nrows : ncols;
    int nnz = 0;  // number of non-zero elements
    for (int i=0; i<nrc; i++) {
        if (sd[i] != 0)
            nnz ++;
    }

    // print nrows (number of rows), ncols (number of columns), and number of non-zero elements
    info = fprintf(mm_file, "%d %d %d\n", nrows, ncols, nnz);
    if (nrows<=0 || ncols<=0) {
        // invalid nrows or ncols
        fclose(mm_file);
        return 4;
    }
    if (info < 0) {
        // file write error
        fclose(mm_file);
        return info;
    }

    // print elements
    for (int i=0; i<nrc; i++) {
        if (sd[i] != 0) {
            if (typeid(Real) == typeid(double))
                info = fprintf(mm_file, "%d %d %.16lg\n", i+1, i+1, sd[i]);
            else  // assume "float"
                info = fprintf(mm_file, "%d %d %g\n", i+1, i+1, sd[i]);
            if (info < 0) {
                // file write error
                fclose(mm_file);
                return info;
            }
        }
    }

    // done writing; now close the file
    if (fclose(mm_file) == EOF) {
        // close file error
        return 2;
    }

    // a successful write, return 0
    return 0;
}


// instantiation by double

template
int matrix_market_diagonal_matrix_write(const int nrows, const int ncols, const double *sd,  // input
                                        const char file_name[],                              // output file name
                                        const char comment[]);                               // optional comment

// instantiation by float

template
int matrix_market_diagonal_matrix_write(const int nrows, const int ncols, const float *sd,   // input
                                        const char file_name[],                              // output file name
                                        const char comment[]);                               // optional comment
