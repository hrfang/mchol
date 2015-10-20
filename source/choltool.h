//----------------------------------------------------------------------------
// Routines for miscellaneous tools for Cholesky or modified Cholesky
// factorization, e.g. permutating rows/columns of a dense matrix.
// This code comes with no guarantee or warranty of any kind.
// coded by hrfang
//----------------------------------------------------------------------------

#ifndef CHOL_TOOL_H
#define CHOL_TOOL_H

// matrices are stored column-wise
// symmetric matrices are stored only with the lower triangular part

// symmetric matrices have the same data structure as the lower triangular matrices
// so the following 2 routines also work for lower triangular matrices

// input symmetric matrix A is of size nrc-by-nrc and stored in sa[]
// output is index j such that A(j+1,j+1) == max{A(i,i):i=1,...,nrc}
template<typename Real>
int get_index_of_max_diagonal_element_of_symmetric_matrix(const int nrc, const Real *sa);

// input symmetric matrix A is of size nrc-by-nrc and stored in sa[]
// output is index j such that |A(j+1,j+1)| == max{|A(i,i)|:i=1,...,nrc}
template<typename Real>
int get_index_of_max_diagonal_magnitude_of_symmetric_matrix(const int nrc, const Real *sa);

// input symmetric matrix A is stored in sa[] as input
// let A2 be the matrix after symmetrically swapping, stored in sa[] as output
// then A2(:,i+1) == A(:,j+1) and A2(:,j+1) == A(:,i+1)
// due to symmetry, A2(i+1,:) == A(j+1,:) and A2(j+1,:) == A(i+1,:)
// let info be the return value
// if info<0, then (-info)th argument has an illegal value
// otherwise, info==0, which means successful exit
// ps. this routine also works for intermediate L in the Cholesky factorization
template<typename Real>
int interchange_row_column_of_symmetric_matrix(const int nrc, Real *sa,
                                               const int i, const int j);

// phase I of the (relaxed) 2-phase strategy:
// the return:
//    let the return value be info
//    if info == -i is negative, then i-th argumement is invalid 
//    otherwise, info is the number of steps processed in phase 1
// arguments:
// 1,2) input symmetric matrix A is of size nrc-by-nrc and stored in sa[]
//    as an in-place update, sa[] also stores the partial factorization on exit
// 3) the factorization is in LDL^T format, with the diagonal matrix D stored in sd[]
// 4) perm is the permutation array that refects pivoting
//    as input, perm[i]==i for i=0,...,nrc-1 is expected, but
//    it can any permutation array reflecting the a priori permutation
// 5) delta>0 is the tolerance
// 6) relax_factor must be in [0,1]; if relax_factor==0, it means not relaxed
//    otherwise, relax_factor is the relax factor (i.e. mu in our 2008 modified Cholesky paper)
template<typename Real>
int phase_one_factorization(const int nrc, Real *sa,
                            Real *sd, int *perm,
                            const Real delta,
                            const Real relax_factor);

#endif  // end of #ifndef CHOL_TOOL_H
