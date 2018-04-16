% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% Solves the L1 regularized least squares problem
%
%  min norm(G*m-d,2)^2+alpha*norm(L*m,1)
%
% using iteratively reweighted least squares.
%
% Inputs:
%       G,d,L,alpha            Problem data.
%       maxiter                (optional) Maximum number of IRLS
%                              iterations. (100)
%       tolx                   (optional) Tolerance on successive
%                              iterates. (1.0e-4)
%       tolr                   (optional) Tolerance below which we consider an
%                              element of L*m to be effectively
%                              zero. (1.0e-4)
%       m0                     (optional) initial solution. (0)
%
% default values are provided for maxiter, tolx, and tolr if they
% are not specified.  You can use [] for any optional parameter to get
% the default value.
%
% Outputs:
%       mreg                   Optimal solution
%       bestobj                Optimal objective value
%       iter                   Number of iterations used
%       objectivevalues        History of objective values.
%
function [mreg, bestobj, iter, objectivevalues] = irlsl1reg(G, d, L, alpha, maxiter, tolx, tolr, m0)
%
% Default for tolr=1.0e-4
%
if ((nargin < 7) || (isempty(tolr)))
    tolr = 1.0e-4;
end
%
% Default for tolx=1.0e-4;
%
if ((nargin < 6) || (isempty(tolx)))
    tolx = 1.0e-4;
end
%
% Default for maxiter=100
%
if ((nargin < 5) || (isempty(maxiter)))
    maxiter = 100;
end
%
% Make space to store the objective values.
%
objectivevalues = zeros(maxiter, 1);
%
% Set tolerances and maximum iterations for lsqr.
%
lsqrtol = 1.0e-6; % This is the regular default for LSQR.
lsqrits = 5 * max(max(size(G))); % Watch for failure of lsqr to converge.
%
% Global variables hold G, L, W, and ALPHA.
%
global IRLSL1REGG;
global IRLSL1REGL;
global IRLSL1REGW;
%
% Initialize the globals.
%
IRLSL1REGG = G;
IRLSL1REGL = L;
IRLSL1REGW = sqrt(alpha/2) * ones(size(L, 1), 1);
%
% Setup up two some useful constants for computing the preconditioner.
% DGTG is the diagonal of G'*G.  LTS is the element-wise square of
% L'.
%
DGTG = zeros(size(G, 2), 1);
for j = 1:size(G, 2)
    DGTG(j) = G(:, j)' * G(:, j);
end
LTS = (L').^2;
%
% Setup the right hand side of the least squares problem.
%
rhs = [d; zeros(size(L, 1), 1)];
%
% Setup the initial solution.
%
if ((nargin < 8) || (isempty(m0)))
    %
    % Start with an initial unweighted solution.
    %
    M = spdiags(sqrt(DGTG+sqrt(alpha/2)*(LTS * IRLSL1REGW)), 0, size(G, 2), size(G, 2));
    [m, flag] = lsqr(@irlsl1regmult, rhs, lsqrtol, lsqrits, M);
    if (flag ~= 0)
        fprintf('lsqr flag=%d\n', flag);
    end
    mreg = m;
    obj = norm(G*m-d)^2 + alpha * norm(L*m, 1);
    bestobj = obj;
    objectivevalues(1) = obj;
else
    m = m0;
    mreg = m;
    obj = norm(G*m-d)^2 + alpha * norm(L*m, 1);
    bestobj = obj;
    objectivevalues(1) = obj;
end
%
% Iterate until maxiter or we converge and return
%
iter = 1;
while (iter < maxiter)
    iter = iter + 1;
    % get get the magnitude of Lm, but don't let any element be less than tolr
    absLm = abs(L*m);
    %
    % A useful diagnostic is to look at the percentage of effectively
    % 0 entries in absLm.  If this is 100%, then we just get least
    % squares and the reweighting isn't effective in getting the
    % 1-norm.  If this is close to 0%, and the weighted least
    % squares problem is well conditioned, we're fine.  However, if
    % this is close to 0% and the weighted least squares problem is
    % poorly conditioned, then tolr should be increased to avoid the
    % ill-conditioning.
    %
    %fprintf('Maximum entry in Lm is %e\n',max(absLm));
    %fprintf('Minimum entry in Lm is %e\n',min(absLm));
    %fprintf('Fraction of effectively 0 entries in Lm is %e\n',...
    %	  sum(absLm<tolr)/length(absLm));
    %
    % Reset anything below tolr to tolr.
    %
    absLm(absLm < tolr) = tolr;
    %
    % build the diagonal weighting matrix for this iteration
    %
    R = 1 ./ absLm;
    IRLSL1REGW = sqrt(alpha/2) * sqrt(R);
    %
    % Store a copy of the previous iterate for convergence test.
    %
    mold = m;
    %
    % get the new iterate
    %
    M = spdiags(sqrt(DGTG+sqrt(alpha/2)*(LTS * IRLSL1REGW)), 0, size(G, 2), size(G, 2));
    [m, flag] = lsqr(@irlsl1regmult, rhs, lsqrtol, lsqrits, M, [], m);
    if (flag ~= 0)
        fprintf('lsqr flag=%d\n', flag);
    end
    %
    % Update the objective.
    %
    obj = norm(G*m-d)^2 + alpha * norm(L*m, 1);
    objectivevalues(iter) = obj;
    if (obj < bestobj)
        bestobj = obj;
        mreg = m;
    end
    %
    % Output the weighted objective.
    %
    %fprintf('Weighted objective is %e\n',obj);
    %
    % Check for convergence.
    %
    if (norm(m-mold) / (1 + norm(mold)) < tolx)
        %
        % Output the number of iterations used.
        %
        fprintf('irlsl1reg used %d iterations\n', iter);
        %
        % Output the optimal objective value.
        %
        fprintf('optimal objective value=%e\n', bestobj);
        %
        % Truncate the list of objective values.
        %
        objectivevalues = objectivevalues(1:iter);
        return
    end
end
%
% We've exceeded the maximum number of iterations.
% Give a warning, if desired, but return best solution.
%
warning('irlslreg1 maximum iterations exceeded.');

