% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% x=irls(A,b,tolr,tolx,p,maxiter)
%
% Uses the iteratively reweight least squares strategy to find an
% approximate L_p solution to Ax=b.
%
% Input Parameters:
%  A       - Matrix of the system of equations.
%  b       - Right hand side of the system of equations.
%  tolr    - Tolerance below which residuals are ignored.
%  tolx    - Stopping tolerance.  Stop when (norm(newx-x)/(1+norm(x)) < tolx)
%  p       - Specifies which p-norm to use (most often, p=1.)
%  maxiter - Limit on number of iterations of IRLS
%
% Output Parameters:
%   x - Approximate L_p solution.
%
function x = irls(A, b, tolr, tolx, p, maxiter)

% Find the size of the matrix A.
[m, n] = size(A);

% Start the first iteration with R=I, and x=A\b (the least
% squares solution)
R = eye(m);
x = A \ b;

% Now loop up to maxiter iterations
iter = 1;
while (iter <= maxiter)
    iter = iter + 1;
    
    % compute the current residual
    r = A * x - b;
    
    % for each row adjust the weighting factor r based on the residual
    for i = 1:m
        if (abs(r(i)) < tolr)
            r(i) = abs(tolr)^(p - 2);
        else
            r(i) = abs(r(i))^(p - 2);
        end
    end
    
    % insert the weighting factors into R
    R = diag(r);
    
    % find the solution to the weighted problem
    newx = (A' * R * A) \ (A' * R * b);
    
    %check for convergence
    if (norm(newx-x) / (1 + norm(newx)) < tolx)
        x = newx;
        return;
    else
        x = newx;
    end
end
warning('irls maximum iterations exceeded.');
