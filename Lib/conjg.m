% Example 10.4 subroutine
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% [xmin,fxmin]=conjg(@f,@g,x0,tol)
%
% Use the method of conjugate gradients (This version uses the 
% method of Fletcher, Reeves, and Powell) to minimize f.
%
%
function [xmin,fxmin]=conjg(func,grad,x0,tol)
global FUN;
global X;
global P;
%
%  Figure out the size of the problem.  
%
n=size(x0,1);
%
%  Initialize x and oldfx.  
%
x=x0;
fx=feval(func,x);
oldfx=fx;
gnew=feval(grad,x);
%
% The main loop.  While the current solution isn't good enough, keep
% trying...  Stop after 500 iterations in the worst case.  Within
% each iteration of the main loop, do up to n iterations of the
% inner loop, but restart after n iterations.  
%
iter=0;
while (((iter == 0) || (abs(fx-oldfx) > tol*(1+abs(fx)))) && (iter < 500))
%
%  Initialize the inner loop.
%
  k=0;
  beta=0;
  p=zeros(n,1);
%
%  The inner loop.
%
  while (((iter == 0) || (abs(fx-oldfx) > tol*(1+abs(fx)))) && (k < n))
%
%  Compute the new search direction.
%
    p=-gnew+beta*p;
%
%  If it isn't a descent direction, then restart the method.
%
    if (gnew'*p >= 0), break; end
%
%  Setup to minimize f(x+alpha*p)
%
    X=x;
    P=p;
    FUN=func;
%
%  Minimize along the line.
%
    [a,u,b]=bracket(@subfun);
    [alphamin,falphamin,newa,newb]=bsearch(@subfun,a,u,b,1.0e-3);
%
%  The new point is at X+alphamin*P.  Update oldx to.
%
    oldx=x;
    oldfx=fx;
    x=X+alphamin*P;
    fx=falphamin;
%
%  Update beta.  
%
    gold=gnew;
    gnew=feval(grad,x);
    beta=(gnew'*gnew)/(gold'*gold);
%
%  Update the iteration counters.
%
    k=k+1;
    iter=iter+1;
%
%  End of the inner loop.
%
  end
%
%  End of the outer loop.
%
end
%
%  Save the results.
%
xmin=x;
fxmin=fx;
%
%  Print out some information
%
fprintf('Done with Conjugate Gradient Method, Iters=%d\n',iter);

