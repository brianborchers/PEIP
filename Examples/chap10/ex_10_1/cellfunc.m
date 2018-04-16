%ray tracing indexing subroutine
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% [cxl,czl] = cellfunc(x1,xv,z1,zv)
%
% find the indices of the cells that are to the upper left of the
% raypath points(x1,z1)
%
% INPUT
%   x1 - a vector of data points along the x axis
%   xv - a vector of the bins along the x axis
%   z1 - a vector of data points along the z axis
%   zv - a vector of the bins along the z axis
%
% OUTPUT
%   cxl - a vector of the index of the greatest element in xv that is less than 
%         x1
%   czl - a vector of the index of the greatest element in zv that is less than 
%         z1

function [cxl,czl] = cellfunc(x1,xv,z1,zv)

cxl=zeros(length(x1),1);
czl=zeros(length(z1),1);
for k=1:length(x1)
    cxl(k)=max([1; find(xv<x1(k))]);
end

for k=1:length(z1)
    czl(k)=max([1; find(zv<z1(k))]);
end
