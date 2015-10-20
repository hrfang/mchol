% we try to minimize the function:
%   f([x1;x2]) = exp(x1)*(x1+x2^2);
% it is defined in myfunc.m, along with its gradient and Hessian
% a (local) minimum is at x=[1;0]

% Newton's method
x0 = [2;3];
fprintf('We minimize a function defined in myfunc.m, starting with the initial x0=[%g;%g].\n', x0(1),x0(2));
[newton_x, newton_fvals, newton_numf, newton_lambdas] = fminunc_newton('myfunc', x0);
fprintf('Newton''s method ends up with x=[%g;%g], f(x)=%g.\n', newton_x(1),newton_x(2), newton_fvals(end));
fprintf('    number of iterations (gradient evaluations, Hessian evaluations): %d\n', length(newton_fvals));
fprintf('    number of function evaluations in line search: %d (total # function evaluations %d)\n', newton_numf, newton_numf+length(newton_fvals));
% newton_fvals(1:end-1)-newton_fvals(2:end) shows the quadratic convergence

% now try the modified Newton method with the following modified Cholesky algorithms.
mchol_methods = {'gmw81', 'gmw1', 'gmw2', 'se90', 'se99', 'se1'};

% note that in the following test, the compiled mex files in "../mex/" are required.
x = {};
fevals = {};
numf = {};
newton_lambdas = {};
for i = 1:length(mchol_methods)
    [x, fvals, numf, lambdas] = fminunc_mnewton('rosenbrock', mchol_methods{i}, x0);
    fprintf('The modified Newton method with %s ends up with x=[%g;%g], f(x)=%g.\n', mchol_methods{i}, x(1),x(2), fvals(end));
    fprintf('    number of iterations (gradient evaluations, Hessian evaluations): %d\n', length(fvals));
    fprintf('    number of function evaluations in line search: %d (total # function evaluations %d)\n', numf, numf+length(fvals));
    mnewton_x{i} = x;
    mnewton_fvals{i} = fvals;
    mnewton_numf{i} = numf;
    mnewton_lambdas{i} = lambdas;
end
% mnewton_fvals{i}(1:end-1)-mnewton_fvals{i}(2:end) shows the quadratic convergence for i = 1:length(mchol_methods)

% for further analysis, if any, the result of Newton's method is stored in newton_x, newton_fvals, newton_numf, newton_lambdas,
% and the results of the modified Newton method are stored in mnewton_x, mnewton_fvals, mnewton_numf, mnewton_lambdas
