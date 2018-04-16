% EM R0 function
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
%  r0=r0(lambda,R,DELTA,H,M,MU,MU0,SIGMA,SIGMA0,D,OMEGA);
%
%  This function computes R0(lambda).  The remaining arguments, R,...,OMEGA
%  describe the physical model including layers, conductivities, and 
%  magnetic permeabilities.
%
function r0=r0(lambda,R,DELTA,H,M,MU,MU0,SIGMA,SIGMA0,D,OMEGA);
%
% Make space for some work arrays and the result.  
%
N=zeros(M,1);
u=zeros(M,1);
Y=zeros(M-1,1);
r0=zeros(length(lambda),1);
%
% Make sure that i is sqrt(-1) in what follows.
%
i=sqrt(-1);
%
% Now, loop through the lambda's.
%
imo=i*MU*OMEGA;
ismo=i*SIGMA.*MU*OMEGA;
for j=1:length(lambda)
%
% First, calculate all of the N and u values.
%
  u=sqrt(lambda(j)^2+ismo);
  N=u./imo;
%
% Now, calculate the Y's.
%
  Y(M-1)=N(M-1)*(N(M)+N(M-1)*tanh(u(M-1)*D(M-1)))/(N(M-1)+N(M)*tanh(u(M-1)*D(M-1)));
  for k=M-2:-1:1
      tanhud=tanh(u(k)*D(k));
      Y(k)=N(k)*(Y(k+1)+N(k)*tanhud)/(N(k)+Y(k+1)*tanhud);
  end
%
% Finally compute r0.
%
  N0=sqrt(lambda(j)^2+i*SIGMA0*MU0*OMEGA)/(i*MU0*OMEGA);
%
% Check the value of r0(lambda)
%
  r0(j)=(N0-Y(1))/(N0+Y(1));
end



