This set of codes (gmw81, gmw1, gmw2, se90, se99, and se1) compute Cholesky
factorizations of real symmetric matrices, modified if necessary to make
them positive definite. They differ in how they compute the modification.

Given an input matrix A, each produces output matrices L, P, and E so that 
            P * (A + E) * P' = L * L'
where P is a permutation matrix, L is lower triangular, and E is the
modification.


======

This folder contains Octave/Matlab scripts to illustrate the modified Newton
methods for unconstrained optimization, by incorporating the modified
Cholesky algorithms into Newton's method.

The 3 key m-files:
------
"fminunc_newton.m" for Newton's method for unconstrained optimization;
"fminunc_mnewton.m" for the modified Newton method for unconstrained
optimization, requiring the mex files compiled in "../mex/";
"linesearch3.m", used by both "fminunc_newton.m" and "fminunc_mnewton.m",
gives a polynomial line search routine for global convergence.

The polynomial line search algorithm implemented in "linesearch3.m" requires
all the search directions being descent. For non-convex problems, Newton's
method may fail when a search direction is not descent in some iteration.
Note that "linesearch3.m" follows much of C. T. Kelley's implementation:
http://www4.ncsu.edu/~ctk/darts/polyline.m
http://www4.ncsu.edu/~ctk/darts/polymod.m

Three example functions are provided: myfunc.m, rosenbrock.m, and griewank.m


%%%%%%

% in Octave/Matlab, try Newton's method for minimizing
%   f([x1;x2]) = exp(x1)*(x1+x2^2);
% this function is defined in myfunc.m, along with its gradient and Hessian.
x0 = [2;3];
[x, fevals, numf, lambdas] = fminunc_newton('myfunc', x0);
x
% on success, we should see that it finds a local minimizer x=[-1;0].
fevals(2:end)-fevals(1:end-1)  % this shows the quadratic convergence.

% now try the modified Newton method with the GMW81 algorithm.
% x0 remains [2;3].
[x, fevals, numf, lambdas] = fminunc_newton('myfunc', 'gmw81', x0);
x
% on success, we should see that it finds a local minimizer x=[-1;0].
fevals(2:end)-fevals(1:end-1)  % this shows the quadratic convergence.
% one can replace 'gmw81' by 'gmw1', 'gmw2', 'se90', 'se99', or 'se1' for
% another modified Cholesky algorithm.

% test_myfunc.m gives a comprehensive test, which shows that in this case,
% Newton's method gives good result, and the modified Newton methods do not
% alter it.

% a sample log file from executing test_myfunc.m in Octave 3.8.1 is provided
% "octave_log/test_myfunc.log"; the corresponding file from Matlab R2015a
% (trial version) is provided "matlab_log/test_myfunc.log".


%%%%%%

% now try to minimize Rosenbrock's function by Newton's method:
%   f([x1;x2]) = (1-x1)^2 + 100*(x2-x1^2)^2;
% it is defined in rosenbrock.m, along with its gradient and Hessian.
% the (global) minimum is at [1;1].
% this function is non-convex and potentially a trouble maker.
x0 = [0;0];
[x, fevals, numf, lambdas] = fminunc_newton('rosenbrock', x0);
x
% on success, we should see that it finds the global minimizer x=[1;1].
fevals(2:end)-fevals(1:end-1)  % this shows the quadratic convergence.
% although this function is non-convex, Newton's method luckily works well as
% the search directions are descent in all directions.

% now try the modified Newton method with the GMW81 algorithm.
% x0 remains [0;0].
[x, fevals, numf, lambdas] = fminunc_mnewton('rosenbrock', 'gmw81', x0);
x
% on success, we should see that it finds the global minimizer x=[1;1].
fevals(2:end)-fevals(1:end-1)  % this shows the quadratic convergence.
% in addition to 'gmw81', one can also try 'gmw1', 'gmw2', 'se90', 'se99',
% or 'se1'.

% test_rosenbrock.m gives a comprehensive test, which shows that in this
% case, Newton's method gives good result, and the modified Newton methods
% do not alter it.

% a sample log file from executing test_rosenbrock.m in Octave 3.8.1 is
% provided "octave_log/test_rosenbrock.log"; the corresponding file from
% Matlab R2015a (trial version) is "matlab_log/test_rosenbrock.log".


%%%%%%

% now try to minimize Griewank's function by Newton's method:
%   f(x) = 1 + (1/4000)*sum_i(x(i)^2) - prod_i(x(i)/sqrt(i))
% it is defined in rosenbrock.m, along with its gradient and Hessian.
% this function has numerous local minima and is notorious for causing
% troubles in optimization.
x0 = [1;1];
[x, fevals, numf, lambdas] = fminunc_newton('griewank', x0);
% Newton's method failed since a non-descent direction is encountered.
% the failure of the assertion qp0<0 (qp0 is the product of the search
% direction and the gradient) says that the search direction is not descent.

% now try the modified Newton method with the GMW81 algorithm.
% x0 remains [1;1].
[x, fevals, numf, lambdas] = fminunc_mnewton('griewank', 'gmw81', x0);
x
fevals(2:end)-fevals(1:end-1)  % this shows the quadratic convergence.
% in addition to 'gmw81', one can also try 'gmw1', 'gmw2', 'se90', 'se99',
% or 'se1'.

% test_griewank.m gives a comprehensive test of the modified Newton methods,
% which shows that in this case, the SE1 algorithm results in a somewhat
% worse local minimum, while the other 5 modified Cholesky algorithms perform
% comparably.

% note that SE90 algorithm is not generally stable, while the other 5
% algorithms usually work. See, e.g., our paper on the modified Cholesky
% algorithms (2008).

% a sample log file from executing test_griewank.m in Octave 3.8.1 is
% provided "octave_log/test_griewank.log"; the corresponding file from
% Matlab R2015a (trial version) is provided "matlab_log/test_griewank.log".
