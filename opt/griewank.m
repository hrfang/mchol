function [f, g, H] = griewank(x, a)
%
%  Usage: [f, g, H] = griewank(x, a);
%
%  The definition of the Griewank function can be found here:
%  http://en.wikipedia.org/wiki/Griewank_function
%  http://mathworld.wolfram.com/GriewankFunction.html
%
%  The Griewank function is defined as:
%  f(x) = 1 + a*sum_i(x(i)^2) - prod_i(x(i)/sqrt(i))
%  If a is not given (nargin<2), the default is a=1/4000.
%
%  Given x, this routine computes the function value f, and
%  optionally the gradient g and Hessian H.
%
%  This function has numerous local minima and is notorious for causing
%  troubles in optimization.

assert(nargin && length(x));

x = x(:);
n = length(x);
if ~exist('a','var')
    a = 1/4000;
end

% compute f(x), the function
cosi = cos(x./sqrt((1:n)'));
pval = prod(cosi);
f = 1 + a*sum(x.*x) - pval;

if nargout>=2
    % compute g(x), the gradient
    tani = tan(x./sqrt((1:n)'));
    a2 = a*2;
    g = zeros(n,1);
    for i = 1:n
        g(i) = a2*x(i) + pval*tani(i)/sqrt(i);
    end
    if nargout>=3
        % compute H(x), the Hessian
        H = zeros(n,n);
        for j = 1:n
            % the diagonal element H(j,j)
            H(j,j) = a2 + pval/j;
            % now the off-diagonal elements
            for i = j+1:n
                H(i,j) = -pval*tani(i)*tani(j) / sqrt(i*j);
                H(j,i) = H(i,j);
            end
        end
    end
end
