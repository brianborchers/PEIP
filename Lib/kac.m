% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% x=kac(A,b,tolx,maxiter)
%
% Implements Kaczmarz's algorithm to solve a system of equations iteratively.
%
% Input Parameters:
%   A - Constraint matrix.
%   b - right hand side.
%   tolx - difference tolerence for successive iterations (stopping criteria).
%   maxiter - maximum iterations (stopping criteria).
%   omegafactor - reduction of omega after each iteration.
%
% Output Parameters:
%      x - solution.
function x = kac(A, b, tolx, maxiter,omegafactor)
%
% Set omegafactor default to 0.95.
%
if (nargin < 5)
  omegafactor=0.95;
end
%
% First, find the size of the matrix.
%
[m, n] = size(A);
%
% Make a copy of A' to speed up some accesses.
%
AT = A';
%
% Setup an initial solution of all zeros.
%
x = zeros(n, 1);
iter = 0;
%
%  Precompute the row norms squared.
%
n2 = zeros(m, 1);
for i = 1:m
    n2(i) = norm(AT(:, i))^2;
end
%
% Initialize omega to 1.0 for the first major iteration.
%
omega=1.0;
%
% The main loop performs iterations of Kaczmarz algorithm until
% maxiters is exceeded or successive iterates differ by less
% than tolx.
%
while (iter <= maxiter)
    % Update the iteration count.
    iter = iter + 1;
    
    %  Start the update cycle with the current solution.
    newx = x;
    
    %  Perform a cycle of m updates.
    for i = 1:m
        newx = newx - omega*((newx' * AT(:, i) - b(i)) / (n2(i))) * AT(:, i);
    end
    
    %  Check for convergence to fixed solution.
    if (norm(newx-x) / (1 + norm(x)) < tolx)
        x = newx;
        return;
    end
    %  Update omega.
    omega=omega*omegafactor;   
    %  Update x for the next major iteration.
    x = newx;
end

% If no convergence
disp('Max iterations exceeded.');
