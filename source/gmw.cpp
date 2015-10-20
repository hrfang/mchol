//----------------------------------------------------------------------------
// Modified Newton methods in the GMW family, including GMW81, GMW-I, GMW-II.
// This code comes with no guarantee or warranty of any kind.
//
// Reference:
// Modified Cholesky Algorithms: A Catalog with New Approaches,
// Haw-ren Fang and Dianne O'Leary,
// Mathematical Programming, Series A, Vol. 115, No. 2, pp. 319--349, 2008
//----------------------------------------------------------------------------

#include <stdlib.h>  // for exit, atof, atoi, atol, srand, rand, etc.
#include <string.h>  // for memcpy, strcmp, strncmp, etc.
#include <stddef.h>  // for NULL pointer, pointer subtraction, etc.
#include <math.h>    // for sqrt, pow, etc.
#include <iostream>  // for cout, cerr, endl, etc. (under namespace std)
#include <limits>    // for std::numeric_limits

#include "choltool.h"
#include "gmw.h"


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
//    delta0 is the modification tolerance if it is greater than 0
//       if delta0==0, then the modification tolerance will be set as the default; delta0<0 is invalid
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
              const Real delta0,
              const int pivot_method,
              const bool is_type1, const bool nondecreasing,
              const bool is_2phase, const Real relax_factor,
              const bool special_last)
{
    // argument verification
    if (nrc <= 0)
        return -1;
    if (sa == NULL)
        return -2;
    if (delta0 < 0.0)
        return -6;
    if (pivot_method<0 || pivot_method>2)
        return -7;
    if (relax_factor<0.0 || relax_factor>1.0)
        return -11;

    const int last_special_steps = special_last ? 2 : 0;  // effective only if (relaxed) 2-phase strategy is applied

    // the diagonal matrix D
    if (sd == NULL)
        sd = new Real[nrc];

    // permutation array for pivoting
    if (perm == NULL) {
        perm = new int[nrc];
        // initialization, no pivoting performed yet
        for (int i=0; i<nrc; i++)
            perm[i] = i;
    }
    // else, assume perm[] already contains a permutation, which can be just the identity, i.e. perm[i]==i

    // diagonal modification
    const bool store_modified = (modified != NULL);
    if (modified == NULL)
        modified = new Real[nrc];

    const Real eps = std::numeric_limits<Real>::epsilon();
    Real delta = delta0;
    if (delta == 0.0) {
        // set delta by the default value
        if (is_type1) {
            // GMW81 or GMW-I algorithm
            delta = eps;
        }
        else {
            // GMW-II algorithm
            Real max_d = 0.0;  // maximum diagonal magnitude of A
            const Real *sa0 = sa;
            for (int i=0; i<nrc; i++) {
                const Real val = (*sa0>=0.0) ? *sa0 : -*sa0;  // abs(A(i+1,i+1)) in Octave/Matlab
                if (max_d < val)
                    max_d = val;
                sa0 += (nrc-i);  // address of the next diagonal element
            }
            // now max_d is max(abs(diag(A))) in Octave/Matlab
            delta = max_d * static_cast<Real>(pow(eps, 2.0/3.0));  // as that used in SE99 algorithm
            if (delta < eps) {
                // in case delta is too small, set it as eps
                delta = eps;
            }
        }
    }

    // phase I, if any
    int steps_processed = 0;
    if (is_2phase) {
        steps_processed = phase_one_factorization(nrc, sa, sd, perm,
                                                  delta, relax_factor);
        for (int i=0; i<steps_processed; i++)
            modified[i] = 0.0;  // no modification in phase I
    }

    // compute beta2 (i.e. beta square), which is max{m1, m2/sqrt(n*n-1), eps} in GMW81, where
    // 1) n is the dimension of n is the dimension of A (input matrix)
    // 2) m1 (eta in the literature) is the maximum diagonal magnitude of A
    // 3) m2 (xi in the literature) is the maximum off-diagonal magnitude of A
    // note that the rationale beta2 >= m1 is to prevent modification when A is sufficiently positive-definite
    // GMW81 is a type I algorithm without using 2-phase or relaxed 2-phase strategy
    // if the (relaxed) 2-phase strategy is applied, then
    // 1) the matrix of concern is now the Schur complement at the beginning of phase II
    // 2) this Schure complement in phase II, if any, needs modified for sure, so beta2 >= m1 is no longer required
    // therefore, we can just use beta2 := max{m2/sqrt(n*n-1), eps}
    // for type II GMW, we replace m2/sqrt(n*n-1) by m2/sqrt(n*n-n)
    // these changes are made with respect to the four objectives; see our paper for details
    Real beta2 = eps;
    Real *sa0 = sa + (2*nrc-steps_processed+1)*steps_processed/2;
    // memory address offset is nrc + (nrc-1) + ... + (nrc-(k-1)) == (2*nrc-k+1)*k/2, where k is steps_processed
    // so sa0 points to A(steps_processed+1,steps_processed+1)
    // update beta2 if is_2phase (2-phase or relaxed 2-phase strategy is applied) and m1>eps
    if (!is_2phase) {
        const Real *sa1 = sa0;
        Real m1 = 0.0;
        for (int i=steps_processed; i<nrc; i++) {
            const Real val = (*sa1>=0.0) ? *sa1 : -*sa1;
            // val is abs(A(i+1,i+1)) in Octave/Matlab, where A(steps_processed:end,steps_processed:end) stores the Schur complement
            if (m1 < val)
                m1 = val;
            sa1 += (nrc-i);  // address of the next diagonal element
        }
        // now m1 is the maximum diagonal magnitude
        if (beta2 < m1)
            beta2 = m1;
    }
    Real m2 = 0.0;
    Real *sa1 = sa0;  // sa0 points to A(i+1,i+1)
    for (int j=steps_processed; j<nrc; j++) {
        sa1 ++;            // skip A(j+1,j+1) in Octave/Matlab
        int ii = nrc-j-1;  // number of off-diagonal elements in (j+1)st column
        while (ii--) {
            const Real val = (*sa1>=0.0) ? *sa1 : -*sa1;
            // val is abs(A(nrc-ii,j+1)) in Octave/Matlab
            if (m2 < val)
                m2 = val;
            sa1 ++;
        }
    }
    // now m2 is the maximum off-diagonal magnitude
    if (is_2phase) {
        const int n2 = nrc-steps_processed-last_special_steps;
        if (n2 > 1) {
            if (is_type1)  // type II
                m2 /= static_cast<Real>(sqrt((Real)(n2*n2-1)));
            else  // type II
                m2 /= static_cast<Real>(sqrt((Real)(n2*n2-n2)));
        }
    }
    else {  // without 2-phase or relaxed 2-phase strategy
        if (nrc > 1) {
            if (is_type1)  // type I
                m2 /= static_cast<Real>(sqrt((Real)(nrc*nrc-1)));
            else  // type II
                m2 /= static_cast<Real>(sqrt((Real)(nrc*nrc-nrc)));
        }
    }
    if (beta2 < m2)
        beta2 = m2;
    Real *arr = new Real[nrc];  // workspace
    int max_steps = nrc;
    if (is_2phase)
        max_steps -= last_special_steps;  // last_special_steps==2 if with the special treatment for the last 1-by-1/2-by-2 matrix
    sa1 = sa0;
    // sa0 points to A(i+1,i+1), the address storing the Schur complement of size (nrc-steps_processed)-by-(nrc-steps_processed)
    for (int i=steps_processed; i<max_steps; i++) {
        // sa1 now points to A(i+1,i+1), the starting element of the Schur complement of size (nrc-i)-by-(nrc-i)
        int idx = i;
        if (pivot_method == 1)
            idx += get_index_of_max_diagonal_element_of_symmetric_matrix(nrc-i, sa1);
        else if (pivot_method == 2)
            idx += get_index_of_max_diagonal_magnitude_of_symmetric_matrix(nrc-i, sa1);
        // idx is the pivot index, indicating the maximum diagonal element (pivot_method==1)
        // or the maximum diagonal magnitude (pivot_method==2)
        if (idx != i) {
            // here i<idx for sure
            // swap (i+1)st and (idx+1)st variables
            interchange_row_column_of_symmetric_matrix(nrc, sa, i, idx);
            // swap perm[i] and perm[idx]
            const int p = perm[i];
            perm[i] = perm[idx];
            perm[idx] = p;
        }
        Real d = *sa1;  // A(i+1,i+1)
        if (is_type1 && d<0.0)
            d = -d;
        // d == |A(i+1,i+1)| if is_type1 == true
        if (d < delta)
            d = delta;
        // now d == max{|A(i+1,i+1)|, delta} if is_type1 == true
        //     d == max{ A(i+1,i+1),  delta} if is_type1 == false
        // compute ||c||_infty, where c is A(i+2:nrc,i+1), i.e. the column under the pivot
        m2 = 0.0;
        Real *sa2 = sa1;  // *sa1 is A(i+1,i+1)
        for (int j=i+1; j<nrc; j++) {
            sa2 ++;
            const Real val = (*sa2>=0.0) ? *sa2 : -*sa2;  // val is |A(j+1,i+1)|
            if (m2 < val)
                m2 = val;
        }
        // now m2 is ||c||_infty, where c is A(i+2:nrc,i+1)
        m2 *= m2;  // so m2 is squared ||c||_infty
        if (m2 > d*beta2)
            d = m2/beta2;
        // now d == max{|A(i+1,i+1)|, delta, (||c||_infty)^2/beta2} if is_type1 == true
        //     d == max{ A(i+1,i+1),  delta, (||c||_infty)^2/beta2} if is_type1 == false
        modified[i] = d-(*sa1);  // *sa1 is A(i+1,i+1), so d == A(i+1,i+1) + modified[i]
        if (nondecreasing && i && modified[i]<modified[i-1]) {
            // increase modified[i] to be modified[i-1]
            d += modified[i-1]-modified[i];
            modified[i] = modified[i-1];
        }
        sd[i] = d;
        // now sd[i] == A(i,i) + modified[i]

        // backup A(i+2:nrc,i+1), the column under pivot
        sa2 = sa1;
        for (int j=i+1; j<nrc; j++) {
            sa2 ++;
            arr[j] = *sa2;  // A(j+1,i+1)
        }

        // compute the (i+1)st column of L
        *sa1 = 1.0;  // A(i+1,i+1) = 1.0
        sa2 = sa1;
        for (int j=i+1; j<nrc; j++) {
            sa2 ++;
            (*sa2) /= d;  // A(j+1,i+1) = A(j+1,i+1) / d
        }

        // update the Schur complement "stored" in A(i+1:nrc,i+1:nrc)
        sa0 = sa1 + (nrc-i);  // pointing to A(i+2,i+2)
        for (int j=i+1; j<nrc; j++) {
            const Real aj = arr[j];
            sa2 = sa1 + (j-i);  // pointer to A(j+1,i+1)
            for (int k=j; k<nrc; k++) {
                // update A(k+1,j+1) = A(k+1,j+1) - arr[j]*A(k+1,i+1),
                // where arr[j] is old A(j+1,i+1) before update and A(k+1,i+1) was already /= A(i+1,i+1)
                *sa0 -= aj*(*sa2);
                // sa0 points to A(k+1,j+1) and sa2 points to A(k+1,i+1)
                // now update sa0 and sa2 for next iteration
                sa0 ++;
                sa2 ++;
            }
        }
        sa1 += (nrc-i);  // this is the next sa1, pointing to A(i+2,i+2)

        // a successful step
        steps_processed ++;
    }
    // now sa1 points to A(steps_processed+1,steps_processed+1), if steps_processed<nrc

    if (steps_processed == nrc-2) {
        // it did get into phase II, with special treatment
        const Real tau = static_cast<Real>(pow(eps, 1.0/3.0));  // eps^(1/3)
        Real &a11 = *sa1, &a21 = *(sa1+1), &a22 = *(sa1+2);
        // denote the last 2-by-2 Schur complement by [ a11 a21 ]
        //                                            [ a21 a22 ]
        Real t1 = a11 + a22;
        Real t2 = a11*a22 - a21*a21;
        // solve x^2 - t1*x + t2 = 0 for eigenvalues of the last 2-by-2 Schur complement lambda1 and lambda2 (lambda1<lambda2)
        t2 = static_cast<Real>(sqrt(t1*t1-4*t2));  // now t2 == lambda2-lambda1,
        Real lambda1 = (t1-t2)/2;                  // and lambda2 is (t1+t2)/2
        t1 = t2 * tau/(1-tau);
        Real s = -lambda1 + (t1>delta ? t1 : delta);
        // we obtain s = -lambda1 + max{(lambda2-lambda1)*tau/(1-tau), delta}

        if (is_type1 && s<-2*lambda1)
            s = -2*lambda1;  // lambda1 < 0.0
        if (s < 0)
            s = 0;  // it would happen only if is_type1==false
        // enforce the nondecreasing strategy, if specified
        if (nondecreasing && modified[nrc-3]>s)
            s = modified[nrc-3];
        modified[nrc-1] = modified[nrc-2] = s;

        // perform modified LDL^T for the last 2-by-2 Schur complement
        sd[nrc-2] = a11 + s;
        a11 = 1.0;
        sd[nrc-1] = a22 + s - a21*(a21/sd[nrc-2]);
        a22 = 1.0;
        a21 /= sd[nrc-2];
        // done
        steps_processed = nrc;
    }
    else if (steps_processed == nrc-1) {
        // it did not get into phase II, with special treatment; it happens only with relaxing
        const Real tau = static_cast<Real>(pow(eps, 1.0/3.0));  // eps^(1/3)
        Real &an = *sa1;  // the last element
        Real t = -an * tau/(1-tau);
        if (t < delta)
            t = delta;
        t -= an;
        // now t = -an + max{-an*tau/(1-tau), delta}
        if (t < 0)
            t = 0;
        if (is_type1 && t<-2*an)
            t = -2*an;
        modified[nrc-1] = t;
        // no need to do anything for enforcing the nondecreasing modification, since modified[i]==0 for i=0,...,nrc-2
        sd[nrc-1] = an + modified[nrc-1];
        an = 1;
        // done
        steps_processed = nrc;
    }
    // else sufficiently positive definite (did not get into phase II) or without special treatment

    if (store_modified) {
        // let modified[] reflect permutation
        for (int i=0; i<nrc; i++)
            arr[perm[i]] = modified[i];
        memcpy(modified, arr, static_cast<size_t>(nrc)*sizeof(Real));  // speedy copy
    }

    // free memory
    delete [] arr;
    if (!store_modified)
        delete [] modified;

    return 0;
}


// the following are 3 wrappers of above routine for the particular members in the GMW family: GMW81, GMW-I, GMW-II

// GMW81 algorithm
template<typename Real>
int mchol_gmw81(const int nrc, Real *sa,
              Real *&sd, int *&perm, Real *modified)
{
    const Real delta = std::numeric_limits<Real>::epsilon();
    const int pivot_method = 2;  // for maximum diagonal magnitude
    const bool is_type1 = true;
    const bool nondecreasing = false;
    const bool is_2phase = false;
    const Real relax_factor = 0.0;  // doesn't matter, since is_2phase == false
    const bool special_last = false;

    // invoke the core routine
    return mchol_gmw(nrc, sa,
                     sd, perm, modified,
                     delta, pivot_method,
                     is_type1, nondecreasing,
                     is_2phase, relax_factor,
                     special_last);
}

// GMW-I algorithm
template<typename Real>
int mchol_gmw1(const int nrc, Real *sa,
              Real *&sd, int *&perm, Real *modified)
{
    const Real delta = std::numeric_limits<Real>::epsilon();
    const int pivot_method = 1;  // for maximum diagonal value
    const bool is_type1 = true;
    const bool nondecreasing = false;
    const bool is_2phase = true;
    const Real relax_factor = 0.75;  // empirically good for the 33 test matrices
    const bool special_last = false;

    // invoke the core routine
    return mchol_gmw(nrc, sa,
                     sd, perm, modified,
                     delta, pivot_method,
                     is_type1, nondecreasing,
                     is_2phase, relax_factor,
                     special_last);
}

// GMW-II algorithm
template<typename Real>
int mchol_gmw2(const int nrc, Real *sa,
              Real *&sd, int *&perm, Real *modified)
{
    // set delta as max_d*eps^(2/3) as that used in the SE99 algorithm, where max_d is the maximum diagonal magnitude of A
    Real max_d = 0.0;  // maximum diagonal magnitude
    const Real *sa0 = sa;
    for (int i=0; i<nrc; i++) {
        const Real val = (*sa0>=0.0) ? *sa0 : -*sa0;  // abs(A(i+1,i+1)) in Octave/Matlab
        if (max_d < val)
            max_d = val;
        sa0 += (nrc-i);  // address of the next diagonal element
    }
    // now max_d is max(abs(diag(A))) in Octave/Matlab
    const Real delta = max_d * static_cast<Real>(pow(std::numeric_limits<Real>::epsilon(), 2.0/3.0));  // as that used in the SE99 algorithm
    const int pivot_method = 1;  // for maximum diagonal value
    const bool is_type1 = false;
    const bool nondecreasing = true;
    const bool is_2phase = true;
    const Real relax_factor = 0.75;  // empirical good for the 33 test matrices
                                     // however, in an experiment of the application to EDMCP, 1.0 works better
    const bool special_last = false;

    // invoke the core routine
    return mchol_gmw(nrc, sa,
                     sd, perm, modified,
                     delta, pivot_method,
                     is_type1, nondecreasing,
                     is_2phase, relax_factor,
                     special_last);
}


// instantiation by double

template
int mchol_gmw(const int nrc, double *sa,
              double *&sd, int *&perm, double *modified,
              const double delta0,
              const int pivot_method,
              bool is_type1, bool nondecreasing,
              bool is_2phase, double relax_factor,
              bool special_last);

template
int mchol_gmw81(const int nrc, double *sa,
              double *&sd, int *&perm, double *modified);

template
int mchol_gmw1(const int nrc, double *sa,
              double *&sd, int *&perm, double *modified);

template
int mchol_gmw2(const int nrc, double *sa,
              double *&sd, int *&perm, double *modified);

// instantiation by float

template
int mchol_gmw(const int nrc, float *sa,
              float *&sd, int *&perm, float *modified,
              const float delta0,
              const int pivot_method,
              bool is_type1, bool nondecreasing,
              bool is_2phase, float relax_factor,
              bool special_last);

template
int mchol_gmw81(const int nrc, float *sa,
              float *&sd, int *&perm, float *modified);

template
int mchol_gmw1(const int nrc, float *sa,
              float *&sd, int *&perm, float *modified);

template
int mchol_gmw2(const int nrc, float *sa,
              float *&sd, int *&perm, float *modified);


//----------------------------------------------------------------------------

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
              const Real delta0,
              const int pivot_method,
              const bool is_type1, const bool nondecreasing,
              const bool is_2phase, const Real relax_factor,
              const bool special_last)
{
    Real *sd = new Real[nrc];
    int info = mchol_gmw(nrc, sa,
                         sd, perm, modified,
                         delta0, pivot_method,
                         is_type1, nondecreasing,
                         is_2phase, relax_factor,
                         special_last);
    if (info < -2) {
        // when info=-i is negative, the i-th argument of the invoked routine mchol_gmw() is invalid
        // the 3rd argument of the invoked mchol_gmw(), sd, is not an argument of the present routine
        info ++;  // so we increase info by 1 to reflect invalid argument of the present routine
    }
    Real *sa0 = sa;
    if (info == 0) {
        // the factorization was successful
        for (int j=0; j<nrc; j++) {
            Real val = static_cast<Real>(sqrt(sd[j]));
            const int nr = nrc-j;
            for (int i=0; i<nr; i++) {
                (*sa0++) *= val;
            }
        }
    }

    // free memory
    delete [] sd;

    return info;
}


// the following are 3 wrappers of the above routine for the particular members in the GMW family: GMW81, GMW-I, GMW-II

// GMW81 algorithm
template<typename Real>
int mchol_gmw81(const int nrc, Real *sa,
                int *&perm, Real *modified)
{
    const Real delta = std::numeric_limits<Real>::epsilon();
    const int pivot_method = 2;  // for maximum diagonal magnitude
    const bool is_type1 = true;
    const bool nondecreasing = false;
    const bool is_2phase = false;
    const Real relax_factor = 0.0;  // doesn't matter, since is_2phase == false
    const bool special_last = false;

    // invoke the core routine
    return mchol_gmw(nrc, sa,
                     perm, modified,
                     delta, pivot_method,
                     is_type1, nondecreasing,
                     is_2phase, relax_factor,
                     special_last);
}

// GMW-I algorithm
template<typename Real>
int mchol_gmw1(const int nrc, Real *sa,
               int *&perm, Real *modified)
{
    const Real delta = std::numeric_limits<Real>::epsilon();
    const int pivot_method = 1;  // for maximum diagonal value
    const bool is_type1 = true;
    const bool nondecreasing = false;
    const bool is_2phase = true;
    const Real relax_factor = 0.75;  // empirical good for the 33 test matrices
    const bool special_last = false;

    // invoke the core routine
    return mchol_gmw(nrc, sa,
                     perm, modified,
                     delta, pivot_method,
                     is_type1, nondecreasing,
                     is_2phase, relax_factor,
                     special_last);
}

// GMW-II algorithm
template<typename Real>
int mchol_gmw2(const int nrc, Real *sa,
               int *&perm, Real *modified)
{
    // set delta as max_d*eps^(2/3) as that used in the SE99 algorithm, where max_d is the maximum diagonal magnitude of A
    Real max_d = 0.0;  // maximum diagonal magnitude
    const Real *sa0 = sa;
    for (int i=0; i<nrc; i++) {
        const Real val = (*sa0>=0.0) ? *sa0 : -*sa0;  // abs(A(i+1,i+1)) in Octave/Matlab
        if (max_d < val)
            max_d = val;
        sa0 += (nrc-i);  // address of the next diagonal element
    }
    // now max_d is max(abs(diag(A))) in Octave/Matlab
    const Real delta = max_d * static_cast<Real>(pow(std::numeric_limits<Real>::epsilon(), 2.0/3.0));  // as that used in the SE99 algorithm
    const int pivot_method = 1;  // for maximum diagonal value
    const bool is_type1 = false;
    const bool nondecreasing = true;
    const bool is_2phase = true;
    const Real relax_factor = 0.75;  // empirical good for the 33 test matrices
                                     // however, in an experiment of the application to EDMCP, 1.0 works better
    const bool special_last = false;

    // invoke the core routine
    return mchol_gmw(nrc, sa,
                     perm, modified,
                     delta, pivot_method,
                     is_type1, nondecreasing,
                     is_2phase, relax_factor,
                     special_last);
}


// instantiation by double

template
int mchol_gmw(const int nrc, double *sa,
              int *&perm, double *modified,
              const double delta0,
              const int pivot_method,
              bool is_type1, bool nondecreasing,
              bool is_2phase, double relax_factor,
              bool special_last);

template
int mchol_gmw81(const int nrc, double *sa,
              int *&perm, double *modified);

template
int mchol_gmw1(const int nrc, double *sa,
              int *&perm, double *modified);

template
int mchol_gmw2(const int nrc, double *sa,
              int *&perm, double *modified);

// instantiation by float

template
int mchol_gmw(const int nrc, float *sa,
              int *&perm, float *modified,
              const float delta0,
              const int pivot_method,
              bool is_type1, bool nondecreasing,
              bool is_2phase, float relax_factor,
              bool special_last);

template
int mchol_gmw81(const int nrc, float *sa,
              int *&perm, float *modified);

template
int mchol_gmw1(const int nrc, float *sa,
              int *&perm, float *modified);

template
int mchol_gmw2(const int nrc, float *sa,
              int *&perm, float *modified);
