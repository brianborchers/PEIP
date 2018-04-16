% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%[mreg,bestobj,iter,objectivevalues]=admml1reg(G,d,L,alpha,maxiter,epsr,lsqrtol,m0)
% Solves the L1 regularized least squares problem
%
%  min norm(G*m-d,2)^2+alpha*norm(L*m,1)
%
% using the Alternating Direction Method of Multipliers (ADMM) method
%
% Inputs:
%       G,d,L,alpha            Problem data.
%       maxiter                (optional) Maximum number of ADMM
%                              iterations. (500)
%       eps                    (optional) Relative error tolerance. (1.0e-5)
%       lsqrtol                (optional) Tolerance for LSQR solves. (1.0e-7)
%       m0                     (optional) initial solution. (0)
%
% default values are provided for maxiter, eps, lsqrtol, and m0 if they
% are not given.  You can use [] for any optional parameter to get
% the default value.
%
% Outputs:
%       mreg                   Optimal solution
%       bestobj                Optimal objective value
%       iter                   Number of iterations used
%       objectivevalues        History of objective values.
%
function [mreg, bestobj, iter, objectivevalues] = ...
    admml1reg(G, d, L, alpha, maxiter, epsr, lsqrtol, m0)
%
% Default for maxiter=50*max(size(G))
%
if ((nargin < 5) || (isempty(maxiter)))
    maxiter = 50 * max(size(G));
end
%
% Default for epsr=1.0e-4.
%
if ((nargin < 6) || (isempty(epsr)))
    epsr = 1.0e-5;
end
%
% Default for lsqrtol=1.0e-6.
%
if ((nargin < 7) || (isempty(lsqrtol)))
    lsqrtol = 1.0e-7;
end
%
% Default for m0 is 0.
%
if ((nargin < 8) || (isempty(m0)))
    m0 = zeros(size(G, 2), 1);
end
%
% Set maximum iterations for lsqr.  This may need to be adjusted
% for problems that are badly conditioned.
%
lsqrits = 5 * max(max(size(G))); % Watch for failure of lsqr to converge.
%
% Parameters for adjusting rho.  If the primal infeasibility is
% more than 10 times smaller than the dual infeasibility or vice
% versa, then we adjust rho by increasing or decreasing it by the
% factor tauinc/taudec.  This helps to balance primal and dual
% feasibility measures.  It's important to not let this change too
% often though, so we don't allow changes on every iteration.
%
%
mu = 10;
tauinc = 2;
taudec = 2;
%
% Make space to store the objective values.
%
objectivevalues = zeros(maxiter, 1);
%
% Global variables hold G, L, W, and ALPHA.
%
global ADMML1REGG;
global ADMML1REGL;
global ADMML1REGRHO;
global ADMML1REGGT;
global ADMML1REGLT;
%
% Initialize the globals.
%
ADMML1REGG = G;
ADMML1REGL = L;
ADMML1REGRHO = 1;
%
% Precompute transposed versions of G and L.
%
ADMML1REGGT = [G', sparse(size(L, 2), size(L, 1))];
ADMML1REGLT = [sparse(size(G, 2), size(G, 1)), L'];
%
% Setup the initial solution.
%
m = m0;
z = zeros(size(L, 1), 1);
u = zeros(size(L, 1), 1);
%
% Compute the initial objective.
%
obj = norm(G*m-d)^2 + alpha * norm(L*m, 1);
objectivevalues(1) = obj;
bestobj = obj;
mreg = m;
%
% Iterate until maxiter or we converge and return
%
iter = 1;
while (iter < maxiter)
    iter = iter + 1;
    %
    % Step 1, update m.
    %
    rhs = [d; sqrt(ADMML1REGRHO/2) * (z - u)];
    %
    % The following alternative code uses direct factorization to
    % solve for m.  It's slow, but can be useful in debugging if lsqr
    % is failing.
    %
    %m=[G; sqrt(ADMML1REGRHO/2)*L]\rhs;
    %
    % Use lsqr() to solve the least squares problem.
    %
    [m, flag] = lsqr(@admml1regmult, rhs, lsqrtol, lsqrits, [], [], m);
    if (flag ~= 0)
        fprintf('lsqr flag=%d\n', flag);
    end
    %
    % Step 2, update z.
    %
    oldz = z;
    z = sthresh(L*m+u, alpha/ADMML1REGRHO);
    %
    % Step 3, update u.
    %
    u = u + (L * m - z);
    %
    % Compute the residuals and report on convergence measures.
    %
    r = L * m - z;
    s = ADMML1REGRHO * L' * (z - oldz);
    normr = norm(r);
    norms = norm(s);
    relnormr = normr / (1 + max(norm(L*m), norm(z)));
    relnorms = norms / (1 + norm(ADMML1REGRHO*(L' * u))); % rho*u=y.
    %
    % Uncomment this line to output the relative errors.
    %
    %fprintf('relnormr=%e, relnorms=%e \n',[relnormr; relnorms]);
    %
    % Update the objective value.
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
    if ((relnormr < epsr) && (relnorms < epsr))
        %
        % We're done.  Clean up the history of objective values and
        % return.
        %
        objectivevalues = objectivevalues(1:iter);
        fprintf('admml1reg used %d iterations\n', iter);
        fprintf('optimal objective value=%e\n', bestobj);
        return
    end
    %
    % periodically, adjust rho if primal or dual infeasibility is out
    % of balance.
    %
    if (mod(iter, 50) == 0)
        %
        % Every 100 iterations, consider adjusting rho.
        %
        if (normr > mu * norms)
            ADMML1REGRHO = tauinc * ADMML1REGRHO;
            u = u / tauinc;
            %
            % For debugging, we can output rho.
            %
            %fprintf('rho=%e\n',ADMML1REGRHO);
        end
        if (norms > mu * normr)
            ADMML1REGRHO = ADMML1REGRHO / tauinc;
            u = u * tauinc;
            %
            % For debugging, we can output rho.
            %
            %fprintf('rho=%e\n',ADMML1REGRHO);
        end
    end
end
%
% We've exceeded the maximum number of iterations.
% Give a warning, if desired, but return best solution.
%
warning('admml1reg maximum iterations exceeded.');

