//----------------------------------------------------------------------------
// Modified Newton methods in the SE family, including SE90, SE99, SE-I.
// This code comes with no guarantee or warranty of any kind.
//
// Reference:
// Modified Cholesky Algorithms: A Catalog with New Approaches,
// Haw-ren Fang and Dianne O'Leary,
// Mathematical Programming, Series A, Vol. 115, No. 2, pp. 319-349, 2008
//----------------------------------------------------------------------------

#include <stdlib.h>  // for exit, atof, atoi, atol, srand, rand, etc.
#include <string.h>  // for memcpy, strcmp, strncmp, etc.
#include <stddef.h>  // for NULL pointer, pointer subtraction, etc.
#include <math.h>    // for sqrt, pow, etc.
#include <iostream>  // for cout, cerr, endl, etc. (under namespace std)
#include <limits>    // for std::numeric_limits

#include "choltool.h"
#include "se.h"

// the modified Cholesky algorithms in the SE family
// given input symmetric matrix A, the output is in the form P*(A+E)*P^T = L*D*L^T, where
//    E is the diagonal modification matrix,
//    P is the permutation matrix, L is unit lower triangular, and D is diagonal,
//    such that A+E is positive definite
// input:
//    A is of size nrc-by-nrc and stored (columnwise) in sa[], only the lower triangular part
// output:
//    the algorithm is in-place; the lower triangular matrix L is stored in sa[];
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
//    gersh_pivot = true for pivoting by the Gershgorin circle theorem
//                = false for no pivoting
//    is_type1 indicates whether the modification is of type 1 or type 2
//    nondecreaseing indicates whether the nondecreasing modification strategy is enforced
//    is_2phase indicates whether the 2-phase strategy is applied or not
//       if 2-phase strategy is applied, then phase 1 is always pivoted by the maximum diagonal element
//    relax_factor is effective only if the 2-phase strategy is applied (i.e. is_2phase==true), in which case
//       if relax_factor==0, it means that the 2-phase strategy is not relaxed
//       otherwise, relax_factor must be in (0,1] as the relaxation factor (i.e. mu in our 2008 modified Cholesky paper)
//    special_last indicates whether the SE special treatment for the last 1-by-1 or 2-by-2 Schur complement is applied
template<typename Real>
int mchol_se(const int nrc, Real *sa,
             Real *&sd, int *&perm, Real *modified,
             const Real delta0,
             const bool gersh_pivot,
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

    Real *gersh = NULL;  // lower end points of Gerschgorin circles, for pivoting
    if (gersh_pivot)
        gersh = new Real[nrc];

    // initialize delta := eps^{2/3} * eta for the SE99 algorithm, or
    // delta := eps^{1/3} * eta for the SE90 algorithm, where eta=max{|a_{ii}|}
    // in case the initial delta < eps, then delta = eps is enforced
    Real eta = 0.0;
    Real *sa0 = sa;
    for (int i=0; i<nrc; i++) {
        const Real val = (*sa0>=0) ? *sa0 : -*sa0;
        if (eta < val)
            eta = val;
        sa0 += (nrc-i);
    }
    const Real eps = std::numeric_limits<Real>::epsilon();
    const Real tau = static_cast<Real>(pow(eps, 1.0/3.0));
    Real delta = delta0;
    if (delta == 0.0) {
        // set delta by the default value
        if (relax_factor>0.0)
            delta = eta*tau*tau;  // used in SE99
        else
            delta = eta*tau;  // used in SE90
        if (delta < eps)
            delta = eps;  // special treatment for both SE90 and SE99
    }

    if (nrc == 1) {  // special trivial case, a 1-by-1 matrix
        Real s = delta-(*sa);  // *sa0 is A(1,1)
        if (is_type1 && s<-2*(*sa))
            s = -2*(*sa);
        if (s<0.0)  // modification must be nonnegative
            s = 0.0;
        // for type I  algorithms, A(1,1)+s = max{delta,|A(1,1)|}
        // for type II algorithms, A(1,1)+s = max{delta, A(1,1) }
        *sd = *sa + s;  // D(1,1) = A(1,1)+s
        *sa = 1.0;      // A(1,1) = 1  (i.e. L(1,1)=1 as output)
        if (store_modified)
            *modified = s;
        else
            delete [] modified;
        return 0;  // success
    }

    // phase I, if any
    int steps_processed = 0;
    if (is_2phase) {
        steps_processed = phase_one_factorization(nrc, sa, sd, perm,
                                                  delta, relax_factor);
        for (int i=0; i<steps_processed; i++)
            modified[i] = 0.0;  // no modification in phase I
    }
    // now steps_processed is the number of steps processed in phase I, which is K_1 in SE's papers (1990,1999)

    if (gersh_pivot) {
        // initialize gersh[] to store the lower end points of Gershgorin circles
        sa0 = sa + (2*nrc-steps_processed+1)*steps_processed/2;
        // memory offset is nrc + (nrc-1) + ... + (nrc-(steps_processed-1)) = (2*nrc-steps_processed+1)*steps_processed/2;
        // sa0 points to A(steps_processed+1,steps_processed+1)
        const int nrc2 = nrc-steps_processed;
        // we consider the nrc2-by-nrc2 Schur complement pointed to by sa0[]
        const Real *sa1 = sa0;
        for (int i=0; i<nrc2; i++) {
            Real g = *sa1;  // A(steps_processed+i+1,steps_processed+i+1), or equivalently
                            // S(i+1,i+1), the (i+1)st diagonal element of the Schur complement
            const Real *sa2 = sa0 + i;
            for (int j=0; j<i; j++) {
                // compute g = g - abs(S(i+1,j+1)), with S the Schur complement
                if (*sa2 > 0)
                    g -= *sa2;
                else
                    g += *sa2;
                sa2 += (nrc2-j-1);
            }
            for (int j=i+1; j<nrc2; j++) {
                // compute g = g - abs(S(j+1,i+1)), with S the Schur complement
                sa2 ++;
                if (*sa2 > 0)
                    g -= *sa2;
                else
                    g += *sa2;
            }
            gersh[steps_processed+i] = g;
            sa1 += (nrc2-i);  // address of the next diagonal element
        }
    }

    int max_steps = nrc;
    if (is_2phase)
        max_steps -= last_special_steps;  // last_special_steps==2 if with the special treatment for 
    Real *arr = new Real[nrc];  // work space
    Real *sa1 = sa0;  // recall that sa0 points to A(steps_processed+1,steps_processed+1)
    for (int i=steps_processed; i<max_steps; i++) {
        // apply the Gershgorin partial pivoting, if gersh_pivot==true
        if (gersh_pivot) {
            int m = i;
            for (int j=i+1; j<nrc; j++) {
                if (gersh[j] > gersh[m])
                    m = j;
            }
            // (m+1)st column has the highest lower Gershgorin end point
            if (m != i) {  // then m>i
                // swap (m+1)st and (i+1)st variables
                interchange_row_column_of_symmetric_matrix(nrc, sa, i, m);
                const Real g = gersh[i];
                gersh[i] = gersh[m];
                gersh[m] = g;
                const int p = perm[i];
                perm[i] = perm[m];
                perm[m] = p;
            }
        }

        // compute the radius of the first Gershgorin circle of the present Schur complement
        Real radius = 0.0;
        // recall that sa1 points to A(i+1,i+1)
        Real *sa2 = sa1;
        const Real aa = *sa2++;  // store the value of A(i+1,i+1) in aa
        for (int j=i+1; j<nrc; j++) {
            // sa2 points to A(i+1,j+1)
            if (*sa2 >= 0)
                radius += *sa2;
            else
                radius -= *sa2;
            sa2 ++;
        }
        // now sa2 points to A(i+2,i+2)

        // compute the modification of the present pivot:
        // s = max{0,-A(i+1,i+1)+max{delta,radius}},               if is_type1 == true
        // s = max{0,-A(i+1,i+1)+max{delta,radius},-2*A(i+1,i+1)}, if is_type1 == false
        Real s = ((radius>=delta) ? radius : delta) - aa;
        if (is_type1 && s<-2*aa)
            s = -2*aa;
        if (s < 0.0)  // modification must be non-negative
            s = 0.0;
        // if nondecreasning==true, s=max{s,modified[i-1]}
        if (nondecreasing && i && modified[i-1]>s)
            s = modified[i-1];
        // update modified[i] and sd[i]
        modified[i] = s;
        sd[i] = aa + modified[i];

        // update the estimated Gershgorin circle lower end points, if the Gershgorin partial pivoting is applied
        if (gersh_pivot) {
            const Real c = 1-radius/sd[i];
            // recall that sa1 points to A(i+1,i+1)
            for (int j=i+1; j<nrc; j++) {
                const Real val = sa1[j-i];  // value of A(j+1,i+1)
                gersh[j] += ((val>=0) ? val : -val)*c;  // increase gersh[j] by abs(val)*c
            }
        }

        // compute (i+1)st column of L
        for (int j=i+1; j<nrc; j++) {
            // recall that sa1 points to A(i+1,i+1)
            arr[j] = sa1[j-i] / sd[i];  // sa1[j-i] is the value of A(j+1,i+1)
        }
        // update for the next Schur complement in the "storage" of A(i+1:n,i+1:n)
        // the pseudo-code is (with ii=i+1):
        // for (int jj=ii+1; jj<=nrc; jj++) {
        //    for (int kk=jj; kk<=nrc; kk++)
        //        A(kk,jj) -= A(kk,ii)*arr[j];  // j=jj-1
        // }
        // the following is a (relatively) efficient implementation
        // recall that sa1 points to A(i+1,i+1) and sa2 points to A(i+2,i+2)
        Real *sa3 = sa2;
        for (int j=i+1; j<nrc; j++) {
            const Real val = arr[j];
            const Real *p = sa1+(j-i);  // p points to A(j+1,i+1)
            for (int k=j; k<nrc; k++) {
                // at the present iteration, sa3 points to A(k+1,j+1) and p points to A(k+1,i+1)
                *sa3++ -= (*p++)*val;
            }
        }

        // update (i+1)st column of L
        *sa1 = 1.0;
        for (int j=i+1; j<nrc; j++)
            sa1[j-i] = arr[j];

        sa1 = sa2;  // now sa1 points to A(i+2,i+2), which is A(i+1,i+1) in the next iteration
        steps_processed ++;  // a successful step
    }
    // if the above for loop was inactive (phase I number of steps >= nrc-2), sa1 points to A(steps_processed+1,steps_processed+1)
    // if the above for loop was active   (phase I number of steps <  nrc-2), sa1 points to A(nrc-1,nrc-1) which remains A(steps_processed+1,steps_processed+1), since steps_processed=nrc-2

    if (steps_processed == nrc-2) {
        // it did get into phase II
        // here sa1 points to A(nrc-1,nrc-1), as steps_processed+1==nrc-1
        // the last 2-by-2 Schur complement is [ A(nrc-1,nrc-1) A(nrc,nrc-1) ], denoted by [ a11 a21 ]
        //                                     [ A(nrc,nrc-1)   A(nrc,nrc)   ]             [ a21 a22 ]
        Real &a11 = *sa1, &a21 = *(sa1+1), &a22 = *(sa1+2);
        Real t1 = a11 + a22;
        Real t2 = a11*a22 - a21*a21;
        // solve x^2 - t1*x + t2=0 for eigenvalues of the last 2-by-2 Schur complement lambda1 and lambda2 (lambda1<lambda2)
        t2 = static_cast<Real>(sqrt(t1*t1-4*t2));  // now t2 == lambda2-lambda1,
        Real lambda1 = (t1-t2)/2;                  // and lambda2 is (t1+t2)/2
        t1 = t2 * tau/(1-tau);
        Real s = -lambda1 + (t1>delta ? t1 : delta);
        // we obtain s = -lambda1 + max((lambda2-lambda1)*tau/(1-tau), delta);

        if (is_type1 && s<-2*lambda1)
            s = -2*lambda1;
        if (s < 0)
            s = 0.0;  // it would happen only if is_type1==false
        // enforce the nondecreasing modification, if specified
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
        // it did not get into phase II; it happens only with relaxing
        // here sa1 points to A(nrc,nrc), as steps_processed+1==nrc
        Real &an = *sa1;  // the last element A(nrc,nrc)
        Real t = -an * tau/(1-tau);
        if (t < delta)
            t = delta;
        t -= an;
        // now t = -an + max{-an*tau/(1-tau), delta}
        if (t < 0.0)
            t = 0.0;
        if (is_type1 && t<-2*an)
            t = -2*an;
        modified[nrc-1] = t;
        // no need to do anything for enforcing the nondecreasing modification, since modified[i]==0 for i=0,...,nrc-2
        sd[nrc-1] = an + modified[nrc-1];
        an = 1.0;
        // done
        steps_processed = nrc;
    }
    // else sufficiently positive definite (did not get into phase II)

    if (store_modified) {
        // let modified[] reflect permutation
        for (int i=0; i<nrc; i++)
            arr[perm[i]] = modified[i];
        memcpy(modified, arr, static_cast<size_t>(nrc)*sizeof(Real));  // speedy copy
    }

    // free memory
    delete [] arr;
    if (gersh_pivot)
        delete [] gersh;
    if (!store_modified)
        delete [] modified;

    return 0;
}


// SE90 algorithm
template<typename Real>
int mchol_se90(const int nrc, Real *sa,
              Real *&sd, int *&perm, Real *modified)
{
    // compute eta, the maximum diagonal magnitude of A
    // in Octave/Matlab, it is eta = max(abs(diag(A))
    Real eta = 0.0;
    const Real *sa0 = sa;
    for (int i=0; i<nrc; i++) {
        const Real val = (*sa0>=0) ? *sa0 : -*sa0;
        if (eta < val)
            eta = val;
        sa0 += (nrc-i);
    }
    // set delta = eta*eps^(1.0/3.0)
    const Real eps = std::numeric_limits<Real>::epsilon();
    const Real tau = static_cast<Real>(pow(eps, 1.0/3.0));
    Real delta = eta*tau;
    if (delta < eps)
        delta = eps;

    // set the other parameters
    const bool gersh_pivot = true;
    const bool is_type1 = false;
    const bool nondecreasing = true;
    const bool is_2phase = true;
    const Real relax_factor = 0.0;  // without relaxation, that causes the instability
    const bool special_last = true;

    // invoke the core routine
    return mchol_se(nrc, sa,
                    sd, perm, modified,
                    delta, gersh_pivot,
                    is_type1, nondecreasing,
                    is_2phase, relax_factor,
                    special_last);
}

// SE99 algorithm
template<typename Real>
int mchol_se99(const int nrc, Real *sa,
              Real *&sd, int *&perm, Real *modified)
{
    // compute eta, the maximum diagonal magnitude of A
    // in Octave/Matlab, it is eta = max(abs(diag(A))
    Real eta = 0.0;
    const Real *sa0 = sa;
    for (int i=0; i<nrc; i++) {
        const Real val = (*sa0>=0) ? *sa0 : -*sa0;
        if (eta < val)
            eta = val;
        sa0 += (nrc-i);
    }
    // set delta = eta*eps^(2.0/3.0)
    const Real eps = std::numeric_limits<Real>::epsilon();
    const Real tau2 = static_cast<Real>(pow(eps, 2.0/3.0));
    Real delta = eta*tau2;
    if (delta < eps)
        delta = eps;

    // set the other parameters
    const bool gersh_pivot = true;
    const bool is_type1 = false;
    const bool nondecreasing = true;
    const bool is_2phase = true;
    const Real relax_factor = static_cast<Real>(0.1);  // >0, with relaxation
    const bool special_last = true;

    // invoke the core routine
    return mchol_se(nrc, sa,
                    sd, perm, modified,
                    delta, gersh_pivot,
                    is_type1, nondecreasing,
                    is_2phase, relax_factor,
                    special_last);
}

// SE-I algorithm
template<typename Real>
int mchol_se1(const int nrc, Real *sa,
              Real *&sd, int *&perm, Real *modified)
{
    // compute eta, the maximum diagonal magnitude of A
    // in Octave/Matlab, it is eta = max(abs(diag(A))
    Real eta = 0.0;
    const Real *sa0 = sa;
    for (int i=0; i<nrc; i++) {
        const Real val = (*sa0>=0) ? *sa0 : -*sa0;
        if (eta < val)
            eta = val;
        sa0 += (nrc-i);
    }
    // set delta = eta*eps^(2.0/3.0)
    const Real eps = std::numeric_limits<Real>::epsilon();
    const Real tau2 = static_cast<Real>(pow(eps, 2.0/3.0));
    Real delta = eta*tau2;
    if (delta < eps)
        delta = eps;

    // set the other parameters
    const bool gersh_pivot = true;
    const bool is_type1 = true;  // the key difference of SE-I from SE90 and SE99
    const bool nondecreasing = false;  // it might be good to set it true
    const bool is_2phase = true;
    const Real relax_factor = static_cast<Real>(0.1);  // >0, with relaxation
    const bool special_last = true;

    // invoke the core routine
    return mchol_se(nrc, sa,
                    sd, perm, modified,
                    delta, gersh_pivot,
                    is_type1, nondecreasing,
                    is_2phase, relax_factor,
                    special_last);
}


// instantiation by double

template
int mchol_se(const int nrc, double *sa,
             double *&sd, int *&perm, double *modified,
             const double delta0,
             const bool gersh_pivot,
             const bool is_type1, const bool nondecreasing,
             const bool is_2phase, const double relax_factor,
             const bool special_last);

template
int mchol_se90(const int nrc, double *sa,
              double *&sd, int *&perm, double *modified);

template
int mchol_se99(const int nrc, double *sa,
              double *&sd, int *&perm, double *modified);

template
int mchol_se1(const int nrc, double *sa,
              double *&sd, int *&perm, double *modified);

// instantiation by float

template
int mchol_se(const int nrc, float *sa,
             float *&sd, int *&perm, float *modified,
             const float delta0,
             const bool gersh_pivot,
             const bool is_type1, const bool nondecreasing,
             const bool is_2phase, const float relax_factor,
             const bool special_last);

template
int mchol_se90(const int nrc, float *sa,
              float *&sd, int *&perm, float *modified);

template
int mchol_se99(const int nrc, float *sa,
              float *&sd, int *&perm, float *modified);

template
int mchol_se1(const int nrc, float *sa,
              float *&sd, int *&perm, float *modified);


//----------------------------------------------------------------------------

// the following routine is a wrapper of the main routine mchol_se(),
// such that the factorization is in the form P*(A+E)*P^T = L*L^T,
// i.e., L is lower triangular but no longer "unit" lower triangular;
// here unit means diagonal elements are all 1s
// the factorization remains in-place and L is stored in sa[]  
// the return:
//    let the return value be info
//    if info == -i is negative, then i-th argumement is invalid
//    otherwise, info == 0
template<typename Real>
int mchol_se(const int nrc, Real *sa,
             int *&perm, Real *modified,
             const Real delta0,
             const bool gersh_pivot,
             const bool is_type1, const bool nondecreasing,
             const bool is_2phase, const Real relax_factor,
             const bool special_last)
{
    Real *sd = new Real[nrc];
    int info = mchol_se(nrc, sa,
                        sd, perm, modified,     
                        delta0, gersh_pivot,   
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

// the following are 3 wrappers of the above routine for the particular members in the SE family: SE90, SE99, SE-I

// SE90 algorithm
template<typename Real>
int mchol_se90(const int nrc, Real *sa,
              int *&perm, Real *modified)
{
    // compute eta, the maximum diagonal magnitude of A
    // in Octave/Matlab, it is eta = max(abs(diag(A))
    Real eta = 0.0;
    const Real *sa0 = sa;
    for (int i=0; i<nrc; i++) {
        const Real val = (*sa0>=0) ? *sa0 : -*sa0;
        if (eta < val)
            eta = val;
        sa0 += (nrc-i);
    }
    // set delta = eta*eps^(1.0/3.0)
    const Real eps = std::numeric_limits<Real>::epsilon();
    const Real tau = static_cast<Real>(pow(eps, 1.0/3.0));
    Real delta = eta*tau;
    if (delta < eps)
        delta = eps;

    // set the other parameters
    const bool gersh_pivot = true;
    const bool is_type1 = false;
    const bool nondecreasing = true;
    const bool is_2phase = true;
    const Real relax_factor = 0.0;  // without relaxation, that causes the instability
    const bool special_last = true;

    // invoke the core routine
    return mchol_se(nrc, sa,
                    perm, modified,
                    delta, gersh_pivot,
                    is_type1, nondecreasing,
                    is_2phase, relax_factor,
                    special_last);
}

// SE99 algorithm
template<typename Real>
int mchol_se99(const int nrc, Real *sa,
              int *&perm, Real *modified)
{
    // compute eta, the maximum diagonal magnitude of A
    // in Octave/Matlab, it is eta = max(abs(diag(A))
    Real eta = 0.0;
    const Real *sa0 = sa;
    for (int i=0; i<nrc; i++) {
        const Real val = (*sa0>=0) ? *sa0 : -*sa0;
        if (eta < val)
            eta = val;
        sa0 += (nrc-i);
    }
    // set delta = eta*eps^(2.0/3.0)
    const Real eps = std::numeric_limits<Real>::epsilon();
    const Real tau2 = static_cast<Real>(pow(eps, 2.0/3.0));
    Real delta = eta*tau2;
    if (delta < eps)
        delta = eps;

    // set the other parameters
    const bool gersh_pivot = true;
    const bool is_type1 = false;
    const bool nondecreasing = true;
    const bool is_2phase = true;
    const Real relax_factor = static_cast<Real>(0.1);  // >0, with relaxation
    const bool special_last = true;

    // invoke the core routine
    return mchol_se(nrc, sa,
                    perm, modified,
                    delta, gersh_pivot,
                    is_type1, nondecreasing,
                    is_2phase, relax_factor,
                    special_last);
}

// SE-I algorithm
template<typename Real>
int mchol_se1(const int nrc, Real *sa,
              int *&perm, Real *modified)
{
    // compute eta, the maximum diagonal magnitude of A
    // in Octave/Matlab, it is eta = max(abs(diag(A))
    Real eta = 0.0;
    const Real *sa0 = sa;
    for (int i=0; i<nrc; i++) {
        const Real val = (*sa0>=0) ? *sa0 : -*sa0;
        if (eta < val)
            eta = val;
        sa0 += (nrc-i);
    }
    // set delta = eta*eps^(2.0/3.0)
    const Real eps = std::numeric_limits<Real>::epsilon();
    const Real tau2 = static_cast<Real>(pow(eps, 2.0/3.0));
    Real delta = eta*tau2;
    if (delta < eps)
        delta = eps;

    // set the other parameters
    const bool gersh_pivot = true;
    const bool is_type1 = true;  // the key difference of SE-I from SE90 and SE99
    const bool nondecreasing = false;  // it might be good to set it true
    const bool is_2phase = true;
    const Real relax_factor = static_cast<Real>(0.1);  // >0, with relaxation
    const bool special_last = true;

    // invoke the core routine
    return mchol_se(nrc, sa,
                    perm, modified,
                    delta, gersh_pivot,
                    is_type1, nondecreasing,
                    is_2phase, relax_factor,
                    special_last);
}


// instantiation by double

template
int mchol_se(const int nrc, double *sa,
             int *&perm, double *modified,
             const double delta0,
             const bool gersh_pivot,
             const bool is_type1, const bool nondecreasing,
             const bool is_2phase, const double relax_factor,
             const bool special_last);

template
int mchol_se90(const int nrc, double *sa,
              int *&perm, double *modified);

template
int mchol_se99(const int nrc, double *sa,
              int *&perm, double *modified);

template
int mchol_se1(const int nrc, double *sa,
              int *&perm, double *modified);

// instantiation by float

template
int mchol_se(const int nrc, float *sa,
             int *&perm, float *modified,
             const float delta0,
             const bool gersh_pivot,
             const bool is_type1, const bool nondecreasing,
             const bool is_2phase, const float relax_factor,
             const bool special_last);

template
int mchol_se90(const int nrc, float *sa,
              int *&perm, float *modified);

template
int mchol_se99(const int nrc, float *sa,
              int *&perm, float *modified);

template
int mchol_se1(const int nrc, float *sa,
              int *&perm, float *modified);
