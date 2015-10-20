//----------------------------------------------------------------------------
// Routines for miscellaneous tools for Cholesky or modified Cholesky
// factorization, e.g. permutating rows/columns of a dense matrix.
// This code comes with no guarantee or warranty of any kind.
// coded by hrfang
//----------------------------------------------------------------------------

#include <stddef.h>  // for NULL pointer, pointer subtraction, etc.
#include <iostream>  // for cout, cerr, endl, etc. (under namespace std)
#include "choltool.h"


// matrices are stored column-wise
// symmetric matrices are stored only with the lower triangular part

// symmetric matrices have the same data structure as the lower triangular matrices
// so the following 2 routines also work for lower triangular matrices


// input symmetric matrix A is of size nrc-by-nrc and stored in sa[]
// output is index j such that A(j+1,j+1) == max{A(i,i):i=1,...,nrc} in Octave/Matlab notation
template<class Real>
int get_index_of_max_diagonal_element_of_symmetric_matrix(const int nrc, const Real *sa)
{
    int j = -1;  // if nrc<1, the return will be -1
    Real maxval = 0.0;
    const Real *sa0 = sa;
    for (int i=0; i<nrc; i++) {
        if (i == 0 || maxval < *sa0) {
            maxval = *sa0;
            j = i;
        }
        sa0 += (nrc-i);
    }

    return j;
}

// input symmetric matrix A is of size nrc-by-nrc and stored in sa[]
// output is index j such that |A(j+1,j+1)| == max{|A(i,i)|:i=1,...,nrc} in Octave/Matlab notation
template<class Real>
int get_index_of_max_diagonal_magnitude_of_symmetric_matrix(const int nrc, const Real *sa)
{
    int j = -1;  // if nrc<1, the return will be -1
    Real maxval = -1.0;
    const Real *sa0 = sa;
    for (int i=0; i<nrc; i++) {
        const Real val = (*sa0>=0.0) ? *sa0 : -*sa0;
        if (maxval < val) {
            maxval = val;
            j = i;
        }
        sa0 += (nrc-i);
    }

    return j;
}

// input symmetric matrix A is stored in sa[] as input
// let A2 be the matrix after symmetrically swapping, stored in sa[] as output
// then A2(:,i+1) == A(:,j+1) and A2(:,j+1) == A(:,i+1)
// due to symmetry, A2(i+1,:) == A(j+1,:) and A2(j+1,:) == A(i+1,:)
// let info be the return value
// if info<0, then (-info)th argument has an illegal value
// otherwise, info==0, which means successful exist
// ps. this routine also works for intermediate L in the Cholesky factorization
template<class Real>
int interchange_row_column_of_symmetric_matrix(const int nrc, Real *sa,
                                               const int i, const int j)
{
    if (nrc <= 0)
        return -1;
    if (sa == NULL)
        return -2;
    if (i<0 || i>=nrc)
        return -3;
    if (j<0 || j>=nrc)
        return -4;
    if (i == j)
        return 0;  // do nothing

    const int ii = (i<j) ? i : j;  // min(i,j)
    const int jj = (i<j) ? j : i;  // max(i,j)
    // we traverse (ii+1)st row with pointer si, and (jj+1)st row with pointer sj
    Real *si = sa + ii;  // address of A(ii+1,1)
    Real *sj = sa + jj;  // address of A(jj+1,1)
    int offset = nrc-1;
    for (int k=0; k<ii; k++) {
        // now *si is A(ii+1,k+1) and *sj is A(jj+1,k+1)
        // swap *si and *sj
        const Real val = *si;
        *si = *sj;
        *sj = val;
        // for the next addresses to swap values
        si += offset;
        sj += offset;
        offset --;
    }

    // now traverse (ii+1)st column with pointer si, and (jj+1)st row with pointer sj
    Real *di = si;  // store the addresss of A(ii+1,ii+1)
    si ++;          // skip A(ii+1,ii+1) to get the address of A(ii+2,ii+1)
    sj += offset;   // skip A(jj+1,ii+1) to get the address of A(jj+1,ii+2)
    offset --;
    for (int k=ii+1; k<jj; k++) {
        // now *si is A(k+1,ii+1) and *sj is A(jj+1,k+1)
        // swap *si and *sj
        const Real val = *si;
        *si = *sj;
        *sj = val;
        // for the next addresses to swap values
        si ++;
        sj += offset;
        offset --;
    }

    // now traverse (ii+1)st column with point si, and (jj+1)st column with pointer sj
    si ++;           // skip A(jj+1,ii+1) to get the address of A(jj+2,ii+1)
    Real *dj = sj;   // store the address of A(jj+1,jj+1)
    sj ++;           // skip A(jj+1,jj+1) to get the address of A(jj+2,jj+1)
    for (int k=jj+1; k<nrc; k++) {
        // now *si is A(k+1,ii+1) and *sj is A(k+1,jj+1)
        // swap *si and *sj
        const Real val = *si;
        *si = *sj;
        *sj = val;
        // for the next addresses to swap values
        si ++;
        sj ++;
    }
    // finally, swap A(ii+1,ii+1) and A(jj+1,jj+1)
    const Real val = *di;
    *di = *dj;
    *dj = val;

    return 0;
}


// instantiation by double

template
int get_index_of_max_diagonal_element_of_symmetric_matrix(const int nrc, const double *sa);

template
int get_index_of_max_diagonal_magnitude_of_symmetric_matrix(const int nrc, const double *sa);

template
int interchange_row_column_of_symmetric_matrix(const int nrc, double *sa,
                                               const int i, const int j);


// instantiation by float

template
int get_index_of_max_diagonal_element_of_symmetric_matrix(const int nrc, const float *sa);

template
int get_index_of_max_diagonal_magnitude_of_symmetric_matrix(const int nrc, const float *sa);

template
int interchange_row_column_of_symmetric_matrix(const int nrc, float *sa,
                                               const int i, const int j);


//----------------------------------------------------------------------------

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
template<class Real>
int phase_one_factorization(const int nrc, Real *sa,
                            Real *sd, int *perm,
                            const Real delta,
                            const Real relax_factor)
{
    // argument verification
    if (nrc <= 0)
        return -1;
    if (sa == NULL)
        return -2;
    if (sd == NULL)
        return -3;
    if (perm == NULL)
        return -4;
    if (delta <= 0.0)
        return -5;
    if (relax_factor<0.0 || relax_factor>1.0)
        return -6;

    // get the threshold
    Real threshold = 0.0;
    if (relax_factor > 0.0) {
        // we first compute the maximum diagonal magnitude of the Schur complement
        Real eta = 0.0;
        const Real *sa0 = sa;
        for (int i=0; i<nrc; i++) {
            // now sa0 points to the address of A(i+1,i+1)
            const Real val = (*sa0>=0.0) ? *sa0 : -*sa0;  // |A(i+1,i+1)|
            if (eta < val)
                eta = val;
            // update sa0 for the address of the next diagonal element
            sa0 += (nrc-i);
        }
        // eta is now the maximum diagonal magnitude of A
        threshold = -relax_factor*eta;  // negative (relaxed)
    }
    else
        threshold = delta;  // positive (not relaxed)

    Real *arr = new Real[nrc];
    Real *sa0 = sa;
    // phase I is always pivoted in this implementation.
    // empirically, pivoting in Phase I improves the performance,
    // in the sense that E (the modification) is reduced
    // note however that it is not required to satisfy the objective 1
    // (i.e., E=0 if A is sufficiently positive definite)
    for (int i=0; i<nrc; i++) {
        // sa0 points to A(i+1,i+1)
        int idx = i + get_index_of_max_diagonal_element_of_symmetric_matrix(nrc-i, sa0);
        if (idx != i) {
            // here i<idx for sure
            // interchange (i+1)st and (idx+1)st variables
            interchange_row_column_of_symmetric_matrix(nrc, sa, i, idx);
            // swap perm[i] and perm[idx]
            const int p = perm[i];
            perm[i] = perm[idx];
            perm[idx] = p;
        }
        Real d = *sa0;  // A(i+1,i+1)
        if (d < delta) {
            delete [] arr;
            return i;  // go to Phase II
        }
        if (relax_factor > 0.0) {
            // relaxed 2-phase strategy
            // we compute the maximum diagonal magnitude of the Schur complement
            Real eta2 = 0.0;
            const Real *sa1 = sa0;  // pointer to A(i+1,i+1)
            for (int j=i; j<nrc; j++) {
                // now sa1 points to A(j+1,j+1)
                const Real val = (*sa1>0) ? *sa1 : -*sa1;  // |A(j+1,j+1)|
                if (eta2 < val)
                    eta2 = val;
                // for the pointer to the next diagonal element
                sa1 += (nrc-j);
            }
            // now eta2 is the maximum diagonal magnitude
            Real threshold2 = -relax_factor*eta2;  // negative
            sa1 = sa0 + (nrc-i);  // pointer to A(i+2,i+2)
            for (int j=i+1; j<nrc; j++) {
                // now sa1 points to A(j+1,j+1)
                if (*sa1 < threshold2) {
                    delete [] arr;
                    return i;  // go to Phase II
                }
                // for the pointer to the next diagonal element
                sa1 += (nrc-j);
            }
        }

        const Real a0 = *sa0;             // A(i+1,i+1)
        const Real *sa1 = sa0 + (nrc-i);  // pointer to A(i+2,i+2)
        Real *sa2 = sa0 + 1;              // pointer to A(i+2,i+1)
        for (int j=i+1; j<nrc; j++) {
            // now sa1 points to A(j+1,j+1), and sa2 points to A(j+1,i+1)
            if (*sa1 - (*sa2)*(*sa2)/a0 < threshold) {
                delete [] arr;
                return i;  // go to Phase II
            }
            // for the pointer to the next diagonal element
            sa1 += (nrc-j);
            // for the pointer to the next off-diagonal element
            sa2 ++;
        }
        // now sa2 points to A(i+2,i+2)

        // update the diagonal element
        sd[i] = a0;  // A(i+1,i+1)
        // compute A(:,i+1), and store it in arr[]
        sa1 = sa0 + 1;  // pointer to A(i+2,i+1)
        for (int j=i+1; j<nrc; j++) {
            // now sa1 points to A(j+1,i+1)
            arr[j] =  *sa1++ / a0;  // store A(j+1,i+1)/A(i+1,i+1) in arr[j]
            // at the same time, sa1 is updated to point to the next off-diagonal element
        }
        // update symmetric A_{22} "stored" in A(i+2:nrc,i+2:nrc)
        // essentially, we update for j=i+1:nrc and k=j:nrc
        //   A(k+1,j+1) = A(k+1,j+1) - A(k+1,i+1) * arr[j]
        // where arr[j] = A(j+1,i+1)/A(i+1,i+1) is pre-computed
        // recall that sa2 points to A(i+2,i+2)
        for (int j=i+1; j<nrc; j++) {
            // sa2 now points to A(j+1,j+1)
            sa1 = sa0 + (j-i);  // address of A(j+1,i+1)
            for (int k=j; k<nrc; k++) {
                // update A(k+1,j+1) = A(k+1,j+1) - A(k+1,i+1)*arr[j]
                *sa2++ -= (*sa1++)*arr[j];
                // at the same time, sa1 and sa2 are updated to point to the next off-diagonal elements
            }
        }

        // update (i+1)st column of A
        *sa0++ = 1.0;  // set A(i+1,i+1)=1.0 and let sa0 points to A(i+2,i+1)
        for (int j=i+1; j<nrc; j++)
            *sa0++ = arr[j];  // set A(j+1,i+1)=arr[j] and let sa0 points to the next element
    }

    // free work space
    delete [] arr;

    return nrc;  // all in phase I, no modifification required
}


// instantiation by float

template
int phase_one_factorization(const int nrc, float *sa,
                            float *sd, int *perm,
                            const float delta,
                            const float relax_factor);

// instantiation by double

template
int phase_one_factorization(const int nrc, double *sa,
                            double *sd, int *perm,
                            const double delta,
                            const double relax_factor);
