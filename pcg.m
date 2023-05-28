function [x, P] = pcg(A, b, tol, iters, B, x0)
%
% Preconditioned Conjugate Gradient method
%
% Input:
%
% - A: The system matrix for the equation Ax = b.
%
% - b: The right-hand side of the equation Ax = b.
%
% - tol: Tolerance for the convergence of the method. If not provided or empty, it defaults to 1e-8.
%
% - iters: The number of iterations to be performed, which is also the number of directions.
%
% - B: preconditioner. If not provided or empty, the standard Conjugate Gradient method is used.
%
% - x0: The initial guess for the solution. If not provided or empty, it defaults to zero.
%
% Output:
%
% - x: The solution to the system Ax = b.
%
% - P: The conjugate directions used in the iterations, normalized to be A-orthonormal.
%
%   This code uses the 'amg' function from the iFEM package developed by L. Chen.
%   Reference:
%   L. Chen. iFEM: an integrated finite element method package in MATLAB. Technical Report, 
%   University of California at Irvine, 2009. 
%   The package is available at: https://github.com/lyc102/ifem

% Default values for optional inputs
if nargin < 3 || isempty(tol), tol = 1e-8; end
if nargin < 6 || isempty(x0), x0 = sparse(length(b), 1); end

n = length(b);
r = b - A * x0; 
x = x0;

% Preconditioner application
if nargin < 5 || isempty(B)  % if B is not supplied, use CG method
    z = r;
else
    if n > 2^10
        option = struct('tol',1e-10,'printlevel',0);
        z = amg(B, r, option);  % for large systems, use algebraic multigrid
    else
        z = B \ r;  % for small systems, solve directly
    end
end

p = z;
rz_prev = dot(r, z);
P = zeros(n, iters + 1);  
P(:, 1) = p;
norm_p = sqrt(dot(A * p, p));
l = zeros(1, iters + 1);  
l(1) = norm_p;

for k = 1:iters - 1
    Ap = A * p;
    alpha = rz_prev / dot(p, Ap);
    x = x + alpha * p;
    r = r - alpha * Ap;
    
    % Preconditioner application
    if nargin < 5 || isempty(B)
        z = r;
    else
        if n > 2^10
            z = amg(B, r, option);
        else
            z = B \ r;
        end
    end

    rz_curr = dot(r, z);
    beta = rz_curr / rz_prev;
    p = z + beta * p;
    P(:, k+1) = p;
    norm_p = sqrt(dot(A * p, p));
    l(k+1) = norm_p;
    rz_prev = rz_curr;

    if sqrt(rz_curr) < tol
        break;
    end
end

l = l(1:k + 1);
P = P(:,1:k + 1); % remove possible redundant zero columns
P = P ./ l;  % normalize columns of P to make it A-orthormal
end
