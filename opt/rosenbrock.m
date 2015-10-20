function [f, g, H] = rosenbrock(x, a, b)
%
%  Usage: [f, g, H] = rosenbrock(x, a, b);
%
%  The definition of the Rosenbrock function can be found here:
%  http://en.wikipedia.org/wiki/Rosenbrock_function
%
%  Given x = [x1; x2], this routine computes the function value, and
%  optionally the gradient and Hessian as:
%    f = (a-x1)^2 + b*(x2-x1^2)^2;
%    g = [ 2*(x1-a) + 4*b*(x1^2-x2)*x1;
%          -2*b*(x1^2-x2) ];
%    H = [ 2+4*b*(3*x1^2-x2), -4*b*x1;
%          -4*b*x1, 2*b ];
%  If a is not given (nargin<2), then the default is a=1.
%  If b is not given (nargin<3), then the default is b=100.
%
%  The (global) minimum is at [a;a^2].
%  This function is non-convex and potentially a trouble maker in
%  optimization.

assert(nargin && length(x)==2);
x1 = x(1);
x2 = x(2);

if ~exist('a','var')
    a = 1;
end
if ~exist('b','var')
    b = 100;
end

% compute f(x), the function
c = x1-a;
d = x1*x1-x2;
f = c*c + b*d*d;

if nargout>=2
    % compute g(x), the gradient
    g = [ 2*c+4*b*d*x1; ...
          -2*b*d ];
    if nargout>=3
        % compute H(x), the Hessian
        H = [ 2+4*b*(3*x1*x1-x2), -4*b*x1; ...
              -4*b*x1, 2*b ];
    end
end
