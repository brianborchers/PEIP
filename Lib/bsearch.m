% Example 10.4 subroutine
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% [xmin,fxmin,newa,newb]=bsearch('f',a,u,b,tol)
%
%
function [xmin,fxmin,newa,newb]=bsearch(fun,a,u,b,tol)
%
%  Useful constants and such. 
%
rooteps=1.0e-8;

%
%  Initialize oldu, oldu2, and oldu3.
%
oldu=u;
oldu2=10*b;
oldu3=20*b;
%
%  Evaluate f at a and b.
%
fa=feval(fun,a);
fb=feval(fun,b);

if (fa < fb)
  x=a;
  fx=fa;
  w=b;
  fw=fb;
else
  x=b;
  fx=fb;
  w=a;
  fw=fa;
end

%
%  evaluate f(u), and update x, w, and v.  
%
fu=feval(fun,u);

if (fu < fx) 
  v=w;
  fv=fw;
  w=x;
  fw=fx;
  x=u;
  fx=fu;
else
  if (fu < fw)
    v=w;
    fv=fw;
    w=u;
    fw=fu;
  end
end

%
%  The main loop.  
%
iter=0;
while ((b-a) > (tol* (1+(a+b)/2)))
%
% update the iteration counter.
%
  iter=iter+1;
  if (iter > 1000)
   disp('bsearch warning: stopping line search after 1000 iterations');
   break; 
  end
%
%  First, compute u, the point the minimizes the quadratic through
%  x, w, and v.  
%

  u=((w^2-v^2)*fx+(v^2-x^2)*fw+(x^2-w^2)*fv)/(2*((w-v)*fx+(v-x)*fw+(x-w)*fv));

%
%  Check to see if u is good enough.  
%

  if ((a < u) && (u < b) && (abs(a-u) > rooteps*(1+abs(a))) &&  (abs(b-u) > rooteps*(1+abs(b))) &&  (abs(u-oldu) < 0.5*abs(oldu2-oldu3)))
%
%  u is ok..
%
%    disp('u is ok ');
%
  else
%
% Need a golden section u.  
%
%    disp('Using the golden section u');

    if ((b-x) > (x-a))
      u=0.3820*(b-x)+x;
    else
      u=0.3820*(x-a)+a;
    end
  end

%
%  Next, compute f(u).
%

  fu=feval(fun,u);

%
%  Update the bracket.
%
  if (x > u) 
    if (fx < fu)
      a=u;
      fa=fu; 
    else
      b=x;
      fb=fx;
    end 
  else
    if (fx < fu)
      b=u;
      fb=fu;
    else
      a=x;
      fa=fx;
    end
  end

%
%  Update x, v, w.
%

  if (fu < fx) 
    v=w;
    fv=fw;
    w=x;
    fw=fx;
    x=u;
    fx=fu;
  else
    if (fu < fw)
      v=w;
      fv=fw;
      w=u;
      fw=fu;
    end
  end

%
%  Update oldu.  
%
  oldu3=oldu2;
  oldu2=oldu;
  oldu=u;

%
% End of the loop.
%
end


xmin=x;
fxmin=fx;
newa=a;
newb=b;





