% Parameter Estimation and Inverse Problems, 2nd edition, 2011 
% by R. Aster, B. Borchers, C. Thurber
%
%
% Construct a two-dimensional, second-order roughening matrix
% for a PSCALE by PSCALE square nodal model (with 
% model wrap-around)
%
%construct roughening matrix
L=zeros(n,n);
k=1;
for j=1:PSCALE
    for i=1:PSCALE;
        mtmp=zeros(PSCALE,PSCALE);
        if i>1; mtmp(i-1,j)=1; end;
        if j>1; mtmp(i,j-1)=1; end;
        if i<PSCALE; mtmp(i+1,j)=1; end;
        if j<PSCALE; mtmp(i,j+1)=1; end;
        L(k,:)=reshape(mtmp,n,1);
        k=k+1;
    end
end

%set values for L diagonal elements, taking into account edges and corners
for i=1:n
    L(i,i)=-sum(L(i,:));
end

