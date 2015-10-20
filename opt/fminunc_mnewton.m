function [x, fvals, numf, lambdas] = fminunc_mnewton(FUNCNAME, MCHOLNAME, x0, tol, maxit, maxnumf)
%
%   The modifed Newton method for minimizing a function.
%
%   Usage: [x, fvals, numf, lambdas] = fminunc_newton(FUNCNAME, MCHOLNAME, x0, tol, maxit, maxnumf)
%
%   We use a modified Cholesky factorization, specified by MCHOLNAME, for
%   computing the modified Newton directions.
%
%   The matlab file FUNCNAME.m defines f(x). FUNCNAME.m takes input x and
%   outputs [f(x),g(x),H(x)], where f(x),g(x),H(x) are the function value,
%   gradient, and Hessian at x.
%
%   The 2nd argument MCHOLNAME specifies the modified Cholesky algorithm.
%   It can be 'gmw81', 'gmw1', 'gmw2', 'se90', 'se99', 'se1'.
%   Note that it is known that 'se90' can result in slow convergence, while
%   the other 5 algorithms are generally stable.
%
%   The 3rd argument x0 is the initial estimate of the minimizer.
%   The 4th argument tol is convergence tolerance (default 1e-6).
%   The 5th agrument maxit is the maximum number of iterations (default 1000).
%   The 6th argument maxnumf is the maximum number of function calls per
%   line search (default 50).
%   numf is the number of function calls.
%
%   The iteration stops when the gradient norm <= tol or when the number of
%   iterations reaches maxit.
%
%   Returned x is the final estimate of the minimizer.
%   fvals is an array of function values f(xi) for i=0,1,... .
%   numf is the number of function calls for line search.
%   lambdas is an array of step lengths at iteration i=0,1,... .
%
%   Example:
%   >> x0 = zeros(2,1);
%   >> [x, fvals] = fminunc_newton('rosenbrock','gmw81',x0)
%   >> fvals(2:end)-fvals(1:end-1)  % this shows the function value change
%   In the example above, f(x) is defined by rosenbrock.m .
%
%   This code comes with no guarantee or warranty of any kind.

path(path, '../mex');  % need compiled gmw81 etc. files in ../mex

assert(nargin>=2);
if ~exist('tol','var') || length(tol)==0
    tol = 1e-6;
end
if ~exist('maxit','var') || length(maxit)==0
    maxit = 1000;
end
if ~exist('maxnumf','var') || length(maxnumf)==0
    maxnumf = 50;
end

x = x0(:);
fvals = [];
numf = 0;
lambdas = [];

% the (modified) Newton iterations
for i = 1:maxit
    [ f, g, H ] = feval(FUNCNAME, x);
    fvals = [ fvals; f ];
    if norm(g) <= tol
        break;
    end
    % perform the modified Choleksy factorization P*(A+E)*P' = L*L'
    [ L, P, E ] = feval(MCHOLNAME, H);
    % solve (A+E)*dx = -g for dx
    % i.e. (P'*L*L'*P)*dx = -g
    dx = -P'*(L'\(L\P*g));
    % line search for step length
    [ x, nf, lam ] = linesearch3(x, f, g, dx, FUNCNAME, maxnumf);
    numf = numf + nf;
    lambdas = [ lambdas; lam ];
    % x = x + dx;  % use this line for no line search
end
