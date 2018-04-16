% Given an image, perform a type-II discrete cosine transform of
% the image, and perform hard threshholding.  Then invert the DCT and
% return the resulting image.
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% [imgout,d]=dctthresh(imgin,threshp)
%
%
% Inputs:
%   imgin           An m by n input image, in gray scale.
%   threshp         fraction of DCT coefficients to zero.
%                   e.g. threshp=0.3 means to 0 out 30% of the coeffs.
%
% Outputs:
%   imgout          The threshholded image.
%   d               The DCT of the with zeroed values
%
function [imgout,d]=dctthresh(imgin,threshp)

% get the size of the image
[m,n]=size(imgin);

% compute its 2D DCT
d=dcost(dcost(imgin.').');

% get the maximum coefficient to zero
dv=reshape(d,m*n,1);
dvs=sort(abs(dv));
thresh=dvs(floor(threshp*m*n));

% zero values below thresh
ind=(abs(d)<=thresh);
d(ind)=0;

% invert the altered DCT
imgout=idcost(idcost(d.').');
