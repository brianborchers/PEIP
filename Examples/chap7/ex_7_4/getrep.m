% Example 7.4
% This routine will create an array of discretized representers
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% Calls the user-supplied kernel function
%	f=getrep(kernfn,xsample,tsample,params);
%
% Inputs
%	kernfn   name of matlab function that calculates the
%	         value of the kernel 
%          This function is called with the
%          following command:
%	         y=kernfn(xx,tt,params)
%	t		     array containing time values
%	xsample  array containing sampling locations
%	tsample  time of sampling
%	params   array of all parameters values
%
% Outputs
%	f		representers

function f=getrep(kernfn,t,xsample,tsample,params);

ns=size(xsample,2);
nt=size(t,2);
weight=t(2)-t(1);
v=params(1);
D=params(2);

f=zeros(ns,nt);

for i=1:nt
	for j=1:ns
		if (t(i) <= tsample) 
			f(j,i)=weight*...
				feval(kernfn,xsample(j),tsample-t(i),params);
		end
	end % for j
end	% for i
