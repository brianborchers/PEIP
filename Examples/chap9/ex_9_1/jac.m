% Example 9.1
% Computes the weighted Jacobian of fun for the model p
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% 
% p is expected to be a two element vector
%   - the first element is transmissivity 
%   - the second element is the storage coefficient
%
function J=jac(p)
% global variables, these are 
% TM    - the times the head was recorded
% SIGMA - the standard deviation for a time
% D     - the distance between the wells
% Q     - the volume of the slug
global TM;
global SIGMA;
global D;
global Q;

% use known formula for the derivatives in the Jacobian
n=length(TM);
J=zeros(n,2);
for i=1:n
  J(i,1)=(-Q*D^2*exp(-D^2*p(1)/(4*p(2)*TM(i)))/...
         (16*pi*p(2)^2*TM(i)^2))/SIGMA(i);
  J(i,2)=(Q/(4*pi*p(2)^2*TM(i)))*...
          ((D^2*p(1))/(4*p(2)*TM(i))-1)*...
          exp(-D^2*p(1)/(4*p(2)*TM(i)))/SIGMA(i);
end
