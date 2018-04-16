% ray tracing subroutine
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% [cxl,czl] = cellfunc(x1,xv,z1,zv)
%
%find the indices of the cells that are to the upper left of the
%raypath points(x1,z1)

function [cxl,czl] = cellfunc(x1,xv,z1,zv)

cxl=ones(length(x1),1);
czl=ones(length(z1),1);
for k=1:length(x1)
    i=find(xv>x1(k),1);
    if ~isempty(i)
    cxl(k)=i-1;
    end
end

for k=1:length(x1)
    i=find(zv>z1(k),1);
    if ~isempty(i)
    czl(k)=i-1;
    end
end
