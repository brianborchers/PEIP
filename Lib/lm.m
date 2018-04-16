% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% [pstar,iter]=lm(func,jac,p0,tol,maxiter)
%
% Use the Levenberg-Marquardt algorithm to minimize
%
%  f(p)=sum(F_i(p)^2)
%
%    func         name of the function F(x)
%    jac          name of the Jacobian function J(x)
%    p0           initial guess
%    tol          stopping tolerance
%    maxiter      maximum number of iterations allowed
%
% Returns
%    pstar        best solution found.
%    iter         Iteration count.
%
function [pstar, iter] = lm(func, jac, p0, tol, maxiter)

%  Initialize p and oldp.
p = p0;
fp = norm(feval(func, p), 2)^2;
oldp = p0 * 2;
oldfp = fp * 2;

% Initialize lambda.
lambda = 0.0001;


% The main loop.  While the current solution isn't good enough, keep
% trying...  Stop after maxiter iterations in the worst case.
%
iter = 0;
while (iter <= maxiter)
    % Compute the Jacobian.
    J = feval(jac, p);
    
    % Compute rhs=-J'*f
    rhs = -J' * feval(func, p);
    
    % Check the termination criteria.
    % all criteria are relative measires
    %
    % the current rhs must be small enough we won't move much AND
    % the change in the norm of f must not have changed to much last step AND
    % the point can't have moved too far last step
    if ((norm(rhs, 2) < sqrt(tol) * (1 + abs(fp))) & ...
            (abs(oldfp-fp) < tol * (1 + abs(fp))) & ...
            (norm(oldp-p, 2) < sqrt(tol) * (1 + norm(p, 2))))
        pstar = p;
        return;
    end
    
    % We use a clever trick here.  The least squares problem
    %
    %  min || [ J              ] s - [ -F ] ||
    %      || [ sqrt(lambda)*I ]     [ 0  ] ||
    %
    % Has the normal equations solution
    %
    %  s=-inv(J'*J+lambda*I)*J'*F
    %
    % which is precisely the LM step.  We can solve this least squares problem
    % more accurately using the QR factorization then by computing
    % inv(J'*J+lambda*I) explicitly.
    %
    myrhs = [-feval(func, p); zeros(length(p), 1)];
    s = [J; sqrt(lambda) * eye(length(p))] \ myrhs;
    
    % compute the new chisq value
    fpnew = norm(feval(func, p+s), 2)^2;
    
    % If the chisq improves then f is improved, make the step, and decrease lambda
    if (fpnew < fp)
        oldp = p;
        oldfp = fp;
        p = p + s;
        fp = fpnew;
        lambda = lambda / 2;
        if (lambda < 10^(-12))
            lambda = 1.0e-12;
        end
    else
        % Didn't improve f, increase lambda, and try again.
        lambda = lambda * 2.5;
        if (lambda > 10^16)
            lambda = 10^16;
        end
    end
    
    %update the iteration count
    iter = iter + 1;
end

%  Return, max iters exceeded.
pstar = p;

