% read the matrix
disp('The matrix:');
A = csvread('first.csv')

% perform the modified Cholesky factorization P*(A+E)*P'=L*L' by the GMW81 algorithm.
[L,P,E] = gmw1(A);
% in addition to gmw81, one can try other algorithms gmw1, gmw2, se90, se99, se1.

% show the modification
disp('The modification by the GMW81 algorithm:');
E

% show residual
res = norm(P*(A+E)*P'-L*L', 'fro');
fprintf('The residual: %g\n', res);

% The matrix stored in first.csv is the first matrix reported in the SE99
% paper that the SE90 algorithm does not work well. One can see this from
% that the modification matrix E by the SE90 algorithm contains much larger
% diagonal elements than the other 5 algorithms.

% Note that the modification matrix E by the 6 algorithms is diagonal and
% nonnegative.
