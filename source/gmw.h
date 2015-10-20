//----------------------------------------------------------------------------
// Modified Newton methods in the GMW family, including GMW81, GMW-I, GMW-II.
// This code comes with no guarantee or warranty of any kind.
//
// Reference:
// Modified Cholesky Algorithms: A Catalog with New Approaches,
// Haw-ren Fang and Dianne O'Leary,
// Mathematical Programming, Series A, Vol. 115, No. 2, pp. 319--349, 2008
//----------------------------------------------------------------------------

#ifndef GMW_H
#define GMW_H

#include <stddef.h>  // for NULL pointer, pointer subtraction, etc.

// the modified Cholesky algorithms in the GMW family
// given input symmetric matrix A, the output is in the form P*(A+E)*P^T = L*D*L^T, where
//    E is the diagonal modification matrix,
//    P is the permutation matrix, L is unit lower triangular, and D is diagonal,
//    such that A+E is positive definite
// input:
//    A is of size nrc-by-nrc and stored (columnwise) in sa[], only the lower triangular part
// output:
//    the algorithm is in-place; the lower triangular matrix L is stored in sa[]
//    the permutation matrix P is represented by the permutation array perm[];
//    if input perm==NULL, then sufficient memory will be allocated and pointed to by perm
//    if input modified!=NULL, the diagonal elements of the modification matrix E
//    are stored in modified[]; otherwise, E will not be stored
// the return:
//    let the return value be info
//    if info==-i is negative, then i-th argumement is invalid
//    otherwise, info==0 and the factorization completes
// algorithm parameters:
//    delta is the modification tolerance if it is greater than 0
//       if delta==0, then the modification tolerance will be set as the default; delta<0 is invalid
//    pivot_method = 0, for no pivoting
//                 = 1, for pivoting by maximum diagonal element
//                 = 2, for pivoting by maximum diagonal magnitude
//    is_type1 indicates whether the modification is of type 1 or type 2
//    nondecreaseing indicates whether the nondecreasing modification strategy is enforced
//    is_2phase indicates whether the 2-phase strategy is applied or not
//       if 2-phase strategy is applied, then phase 1 is always pivoted by the maximum diagonal element
//    relax_factor is effective only if the 2-phase strategy is applied (i.e. is_2phase==true), in which case
//       if relax_factor==0, it means that the 2-phase strategy is not relaxed
//       otherwise, relax_factor must be in (0,1] as the relaxation factor (i.e. mu in our 2008 modified Cholesky paper)
//    special_last indicates whether the SE special treatment for the last 1-by-1 or 2-by-2 Schur complement is applied
template<typename Real>
int mchol_gmw(const int nrc, Real *sa,
              Real *&sd, int *&perm, Real *modified,
              const Real delta,
              const int pivot_method,
              const bool is_type1, const bool nondecreasing,
              const bool is_2phase, const Real relax_factor,
              const bool special_last);

// the following are 3 wrappers of the above routine for the particular members in the GMW family: GMW81, GMW-I, GMW-II

// GMW81 algorithm
template<typename Real>
int mchol_gmw81(const int nrc, Real *sa,
              Real *&sd, int *&perm, Real *modified);

// GMW-I algorithm
template<typename Real>
int mchol_gmw1(const int nrc, Real *sa,
              Real *&sd, int *&perm, Real *modified);

// GMW-II algorithm
template<typename Real>
int mchol_gmw2(const int nrc, Real *sa,
              Real *&sd, int *&perm, Real *modified);


// the following routine is a wrapper of the main routine mchol_gmw(),
// such that the factorization is in the form P*(A+E)*P^T = L*L^T,
// i.e., L is lower triangular but no longer "unit" lower triangular;
// here unit means diagonal elements are all 1s
// the factorization remains in-place and L is stored in sa[]
// the return:
//    let the return value be info
//    if info == -i is negative, then i-th argumement is invalid
//    otherwise, info == 0
template<typename Real>
int mchol_gmw(const int nrc, Real *sa,
              int *&perm, Real *modified,
              const Real delta,
              const int pivot_method,
              const bool is_type1, const bool nondecreasing,
              const bool is_2phase, const Real relax_factor,
              const bool special_last);

// the following are 3 wrappers of the above routine for the particular members in the GMW family: GMW81, GMW-I, GMW-II

// GMW81 algorithm
template<typename Real>
int mchol_gmw81(const int nrc, Real *sa,
              int *&perm, Real *modified);

// GMW-I algorithm
template<typename Real>
int mchol_gmw1(const int nrc, Real *sa,
              int *&perm, Real *modified);

// GMW-II algorithm
template<typename Real>
int mchol_gmw2(const int nrc, Real *sa,
              int *&perm, Real *modified);

#endif  // end of #ifndef GMW_H
