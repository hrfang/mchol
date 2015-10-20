function [f, g, H] = myfunc(x)
%
%  Usage: [f, g, H] = myfunc(x);
%
%  Given x = [x1; x2], this routine computes the function value at x as:
%    f = exp(x1)*(x1+x2^2);
%
%  it also optionally computes the gradient (nargout>=2) and Hessian
%  (nargout>=3) as:
%    g = exp(x1)*[ x1+x2^2+1;
%                  2*x2 ];
%    H = exp(x1)*[ x1+x2^2+2, 2*x2;
%                  2*x2, 2 ];

assert(nargin && length(x)==2);
x1 = x(1);
x2 = x(2);

a = exp(x1);
b = x1+x2*x2;
f = a*b;
if nargout>=2
    g = a*[ b+1; ...
            2*x2 ];
    if nargout>=3
        H = a*[ b+2, 2*x2; ...
                2*x2, 2 ];
    end
end
