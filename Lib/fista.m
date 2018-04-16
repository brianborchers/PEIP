% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% [mreg,bestobj,iter]=ista(G,d,alpha,Lip,maxiter,tolx,m0)
%
% Minimizes
%   norm(G*m-d,2)^2+alpha*norm(m,1)
%
% by fast iterative soft threshholding.
%
% Inputs:
%       G,d,alpha            Problem data.
%       Lip                  Lipschitz constant. (must have L > 2*norm(G'*G,2)
%       maxiter              (optional) Maximum iterations.
%       tolx                 (optional) Convergence tolerance. (1.0e-4)
%       m                    (optional) initial solution.
%
% Outputs:
%       mreg                 Regularized solution
%       bestobj              Objective value of mreg
%       iter                 Number of iterations used.
%
% Notes:
%   Lip must be greater than 2*norm(G'*G).  Since
%   2*norm(G'*G)<2*norm(G)*norm(G'), a safe choice is
%   Lip=2.05*normest(G)*normest(G')
%
%   The default value of maxiter is 20 times the maximum dimension
%   of G.
%
%   The termination criterion is
%      norm(m-oldm)/(1+norm(oldm)) < tolx
%
%   where the default value of tolx is 1.0e-4.
%
function [mreg, bestobj, iter] = fista(G, d, alpha, Lip, maxiter, tolx, m0)
%
% Default values for parameters.
%
if ((nargin < 6) || (isempty(tolx)))
    tolx = 1.0e-4;
end
if ((nargin < 5) || (isempty(maxiter)))
    maxiter = 50 * max(max(size(G)));
end
%
% Initialize the solution.
%
if ((nargin < 7) || (isempty(m0)))
    m = zeros(size(G, 2), 1);
else
    m = m0;
end
%
% Initialize the objective value.
%
y = m;
obj = norm(G*m-d)^2 + alpha * norm(m, 1);
bestobj = obj;
mreg = m;
%
% Initialize the extrapolation parameter theta.
%
theta = 1;
%
% Main loop.
%
iter = 1;
while (iter < maxiter)
    iter = iter + 1;
    %
    % Store a copy of m for checking convergence.m
    %
    mold = m;
    %
    % Compute the gradient of norm(G*y-d)^2
    %
    delf = 2 * G' * (G * y - d);
    %
    % Compute m by soft threshholding.
    %
    m = sthresh(y-(1 / Lip)*delf, alpha/Lip);
    %
    % Update theta.
    %
    thetaold = theta;
    theta = (1 + sqrt(1+4*theta^2)) / 2;
    %
    % Update y.
    %
    y = m + ((thetaold - 1) / theta) * (m - mold);
    %
    % Update the objective.
    %
    obj = norm(G*m-d)^2 + alpha * norm(m, 1);
    %
    % Keep track of the best solution.
    %
    if (obj < bestobj)
        bestobj = obj;
        mreg = m;
    end
    %
    % Check for convergence.
    %
    if (norm(m-mold) / (1 + norm(mold)) < tolx)
        fprintf('FISTA used %d iterations\n', iter);
        return
    end
    %
    % Print out the objective value.
    %
    if (mod(iter, 500) == 0)
        fprintf('Iter=%d, Conv=%e, Obj=%e\n', [iter; norm(m-mold) / (1 + ...
            norm(mold)); obj]);
    end
end
%
% We've exceeded the maximum number of iterations.
% Give a warning, if desired, but return best solution.
%
warning('FISTA maximum iterations exceeded.');
