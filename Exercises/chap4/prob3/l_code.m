% Parameter Estimation and Inverse Problems, 2nd edition, 2011 
% by R. Aster, B. Borchers, C. Thurber
%
%L matrix generating code for Chapter 4, Problem 3
L=zeros(14*14,256);
k=1;
for i=2:15,
for j=2:15,
M=zeros(16,16);
M(i,j)=-4;
M(i,j+1)=1;
M(i,j-1)=1;
M(i+1,j)=1;
M(i-1,j)=1;
L(k,:)=reshape(M,256,1)';
k=k+1;
end
end
