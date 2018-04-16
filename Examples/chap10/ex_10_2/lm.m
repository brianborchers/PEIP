% Levenberg-Marquardt function
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% [xstar,iter]=lm(func,jac,x0,tol,maxiter)
%
% Use the Levenberg-Marquardt algorithm to minimize 
%
%  f(x)=sum(F_i(x)^2)
%
%    func         name of the function F(x)
%    jac          name of the Jacobian function J(x)
%    x0           initial guess
%    tol          stopping tolerance
%    maxiter      maximum number of iterations allowed
%
% Returns
%    xstar        best solution found. 
%    iter         Iteration count.
%
function [pstar,iter]=lm(func,jac,p0,tol,maxiter)
%
%  Initialize p and oldp.  
%
p=p0;
fp=norm(feval(func,p),2)^2;
oldp=p0*2;
oldfp=fp*2;
%
% Initialize lambda.
%
lambda=0.0001;
%
% The main loop.  While the current solution isn't good enough, keep
% trying...  Stop after maxiter iterations in the worst case.
%
  iter=0;
  while (iter <= maxiter)
%
% Compute the Jacobian.
%
    J=feval(jac,p);
%
% Compute rhs=-J'*f
%
%    rhs=-J'*feval(func,p);
%
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
    rhs=-J'*feval(func,p);
    myrhs=[-feval(func,p); zeros(length(p),1)];
    s=[J; sqrt(lambda)*eye(length(p))]\myrhs;
%
% Check the termination criteria.
%
    if ((norm(rhs,2)< sqrt(tol)*(1+abs(fp))) & ...
        (abs(oldfp-fp)<tol*(1+abs(fp))) & ...
        (norm(oldp-p,2)<sqrt(tol)*(1+norm(p,2))))
      pstar=p;
      return;
    end
%
% See whether this improves chisq or not.
%
    fpnew=norm(feval(func,p+s),2)^2;
%
% If this improves f, then make the step, and decrease lambda and make
% the step.
%
    if (fpnew < fp)
      oldp=p;
      oldfp=fp;
      p=p+s;
      fp=fpnew;
      lambda=lambda/2;
      if (lambda <10^(-12))
        lambda=1.0e-12;
      end
      iter=iter+1;
    else
%
% Didn't improve f, increase lambda, and try again.
%
      lambda=lambda*2.5;
      if (lambda >10^16)
        lambda=10^16;
      end
      iter=iter+1;
    end
%
%  end of the loop.
%
  end
%
%  Return, max iters exceeded.
%
pstar=p;




