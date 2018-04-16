% line search function
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
%
%
% [xmin,fxmin,newa,newb]=gssearch('f',a,u,b,tol)
%
%
%
%
function [xmin,fxmin,newa,newb]=gssearch(fun,a,u,b,tol)
%
%  Evaluate f at the end points of the interval, and initialize
%  xmin.
%
newa=a;
fa=feval(fun,newa);
newb=b;
fb=feval(fun,newb);
if (fa < fb)
  xmin=newa;
  fxmin=fa;
else
  xmin=newb;
  fxmin=fb;
end;

%
% Calculate the locations of the two points inside the interval.
%
u=newa+.3820*(newb-newa);
v=newa+0.6180*(newb-newa);

%
%  Calculate f(u)
%
fu=feval(fun,u);

%
% Update xmin and fxmin.  
%
if (fu < fxmin) 
  xmin=u;
  fxmin=fu;
end;

%
%  Calculate f(v)
%

fv=feval(fun,v);

%
%  Update fxmin and xmin.
%
if (fv < fxmin) 
  xmin=v;
  fxmin=fv;
end;

%
%  Now, the main loop.  
%
while ((newb-newa) > (tol* (1+(newa+newb)/2)))

%
% Decide whether to move to [a,v] or [u,b]
%
  if (fu < fv)
%
%  New interval [a,v]
%
    newb=v;
    fb=fv;
%
%  Find the new v and u and evaluate the functions.
%
    v=u;
    fv=fu;
    u=newa+0.3820*(newb-newa);
    fu=feval(fun,u);
%
%  Update xmin and fxmin.
%
    if (fu < fxmin) 
      xmin=u;
      fxmin=fu;
    end;
  else
%
%  New interval [u,b]
%
    newa=u;
    fa=fu;
%
%  Find the new v and u and evaluate the functions.
%
    u=v;
    fu=fv;
    v=newa+0.6180*(newb-newa);
    fv=feval(fun,v);
%
%  Update xmin and fxmin.
%
    if (fv < fxmin) 
      xmin=v;
      fxmin=fv;
    end;
%
% End the while loop.
%
  end;
%
%  End the whole function.
%
end;
