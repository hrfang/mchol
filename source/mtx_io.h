//----------------------------------------------------------------------------
// Routines for reading/writing dense matrix files in matrix-market format.
// This code comes with no guarantee or warranty of any kind.
// coded by hrfang
//----------------------------------------------------------------------------

#ifndef MTX_IO_H
#define MTX_IO_H

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
template<typename Real>
int matrix_market_matrix_read(const char file_name[],                         // input
                              int &nrows, int &ncols, Real *&sa, char &symm,  // output
                              const bool convert=false);                      // optional parameter

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
template<typename Real>
int matrix_market_matrix_write(const int nrows, const int ncols, const Real *sa, const char symm,  // input
                               const char file_name[],                                             // output file name
                               const char comment[]=NULL);                                         // optional comment

// matrix_market_diagonal_matrix_write() writes a diagonal matrix D to file_name[] in matrix-market sparse matrix format
//    where the diagonal of D is stored in sd[], and
//    nrows, ncols are the number of rows, the number of columns, respectively, and
// the return value:
//    0: a successful write
//    1: file open error
//    2: file close error
//    3: invalid nrc (i.e. not positive)
//    a negative integer: a number returned by fprintf() which signals an error
template<typename Real>
int matrix_market_diagonal_matrix_write(const int nrows, const int ncols, const Real *sd,  // input
                                        const char file_name[],                            // output file name
                                        const char comment[]=NULL);                        // optional comment


#endif  // end of #ifdef MTX_IO_H
