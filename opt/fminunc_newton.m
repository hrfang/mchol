function [x, fvals, numf, lambdas] = fminunc_newton(FUNCNAME, x0, tol, maxit, maxnumf)
%
%   Newton's method for minimizing a function.
%
%   Usage: [x, fvals, numf, lambdas] = fminunc_newton(FUNCNAME, x0, tol, maxit);
%
%   The matlab file FUNCNAME.m defines f(x). FUNCNAME.m takes input x and
%   outputs [f(x),g(x),H(x)], where f(x),g(x),H(x) are the function value,
%   gradient, and Hessian at x.
%
%   The 2nd argument x0 is the initial estimate of the minimizer.
%   The 3rd argument tol is convergence tolerance (default 1e-6).
%   The 4th agrument maxit is the maximum number of iterations (default 500).
%   The 5th argument maxnumf is the maximum number of function calls per
%   line search (default 50).
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
%   >> [x, fvals] = fminunc_newton('myfunc', x0)
%   >> fvals(2:end)-fvals(1:end-1)  % this shows the function value change
%   In the example above, f(x) is defined by myfunc.m .
%
%   This code comes with no guarantee or warranty of any kind.


assert(nargin>=2);

if ~exist('tol','var') || length(tol)==0
    tol = 1e-6;
end
if ~exist('maxit','var') || length(maxit)==0
    maxit = 500;
end
if ~exist('maxnumf','var') || length(maxnumf)==0
    maxnumf = 50;
end

x = x0(:);
fvals = [];
numf = 0;
lambdas = [];

% the Newton iterations
for i = 1:maxit
    [ f, g, H ] = feval(FUNCNAME, x);
    fvals = [ fvals; f ];
    if norm(g) <= tol
        break;
    end
    dx = -(H\g);
    [ x, nf, lam ] = linesearch3(x, f, g, dx, FUNCNAME, 20);
    numf = numf + nf;
    lambdas = [ lambdas; lam ];
    % x = x + dx;  % use this line for no line search
end
