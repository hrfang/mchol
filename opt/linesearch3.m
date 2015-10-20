function [xt, numf, lamc] = linesearch3(xc, fc, gc, d, FUNCNAME, maxnumf)
%
%   Polynomial line search.
%
%   Usage: [xt, numf, lambda] = linesearch3(xc, fc, gc, d, FUNCNAME, maxnumf);
%
%   On success, the sufficient decrease condition (a.k.a. Armoji rule) is
%   satisfied.
%
%   The line search method is described in both:
%   [1] "Iterative Methods for Optimization" by C.T. Kelley (chapter 3.2.1)
%   [2] "Numerical Optimization" by J. Nocedal and S.J. Wright (chapter 3.5)
%   The code follows much of C. T. Kelley's implementation:
%   http://www4.ncsu.edu/~ctk/darts/polyline.m
%   http://www4.ncsu.edu/~ctk/darts/polymod.m
%
%   In short, if the Armoji rule (i.e. sufficient decrease condition) does
%   not hold with the full step length, then the minimizer of a quadratic
%   model is computed in the first iteration. If it does not lead to a step
%   length satisfying the Armoji rule, then the iteration continues with
%   cubic models, until a step length is accepted or the maximum number of
%   function evaluations is reached.
%
%   Input:
%   xc = current point
%   fc = current function value
%   gc = current gradient value
%   d  = search direction
%   FUNCNAME = objective function
%              the function will be called by fout = feval(FUNCNAME, x),
%              where the output fout is a real-valued scalar
%   maxnumf  = maximum number of step length reductions
%
%   Output:
%   xt   = new point
%   numf = number of calls to FUNCNAME, if line search succeeds; or
%          maxnumf+1, if maximum number of step length reductions is reached
%   lambda = step length (i.e. the new point is xt = xc + lambda*d)
%
%   This code comes with no guarantee or warranty of any kind.


assert(maxnumf>=1);

% line search parameters
alp       = 1e-4;
beta_low  = 0.1;
beta_high = 0.5;

% first step
lamc = 1;  % initial full step length is typical for Newton-type methods (Newton, quasi-Newton, modified Newton)
xt   = xc + lamc*d;
ft   = feval(FUNCNAME, xt);
numf = 1;  % number of function calls

% set up parameters for the first (quadratic) model/approximation
q0   = fc;
qp0  = gc'*d;     % must be negative (i.e. d is a descent direction), or it doesn't work
assert(qp0 < 0);  % make sure it is a descent direction
qc   = ft;

fgoal = fc + alp*lamc*qp0;
while ft > fgoal && numf <= maxnumf
    % ft <= fgoal is the sufficient decrease condition, known as the Armoji rule
    if numf == 1  % 1st iteration, quadratic
        lam = quadratic_step(q0, qp0, lamc, qc);
    else  % sine the 2nd iteration, cubic
        lam = cubic_step(q0, qp0, lamc, qc, lamm, qm);
    end
    % if lam is not in [lamc*beta_low,lamc*beta_high], project it to this region
    if lam < lamc*beta_low
        lam = lamc*beta_low;
    elseif lam > lamc*beta_high
        lam = lamc*beta_high;
    end
    % lam now is the newly estimated step length
    qm = qc;
    lamm = lamc;
    lamc = lam;
    xt = xc + lamc*d;
    ft = feval(FUNCNAME, xt);
    numf = numf+1;
    qc = ft;
    fgoal = fc + alp*lamc*qp0;
end


function lambda = quadratic_step(q0, qp0, lamc, qc)
%
%   Quadratic polynomial line search step.
%
%   Usage: lambda = quadratic_step(q0, qp0, qc, beta_low, beta_high)
%
%   This routine finds the minimizer of the quadratic polynomial q(lam) that
%   satisfies:
%
%   q(0)=q0, q'(0)=qp0, q(lamc)=qc
%
%   Note that this search step is performed due the the failure of the
%   sufficient decrease condition, which guarantees that the quadratic
%   polynomial is convex provided that qp0<0 (i.e. the search direction is
%   descent) and lamc>0.

den = 2*(qc-q0-qp0*lamc);
lambda = -qp0*lamc*lamc / den;


function lambda = cubic_step(q0, qp0, lamc, qc, lamm, qm)
%
%   Cubic polynomial line search step.
%
%   Usage: lambda = cubic_step(q0, qp0, qc, lamm, qm)
%
%   This routine finds the minimizer of the cubic polynomial q(lam) for
%   lam>0, with q satisfying:
%
%   q(0)=q0, q'(0)=qp0, q(lamc)=qc, q(lamm)=qm

A = [ lamc^2, lamc^3; lamm^2, lamm^3 ];
b = [ qc-q0-qp0*lamc; qm-q0-qp0*lamm ];
c = A\b;
lambda = (-c(1) + sqrt(c(1)*c(1)-3*c(2)*qp0)) / (3*c(2));
