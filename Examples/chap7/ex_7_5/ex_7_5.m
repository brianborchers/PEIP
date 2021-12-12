% Example 7.5.
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% Reset the RNG's.
%
rand('seed',0);
randn('seed',0);
%
% read in the raw image and make it floating point.
%
img=double(imread('image.png'));
%
% Get size parameters.
%
[rows,cols]=size(img);
%
% Plot it.
%

disp('Note: This example may take several minutes to run on most computers.')

figure(1);
clf
imagesc(img);
colormap gray;
bookfonts
axis off

disp('Displaying the original image (fig. 1)')
print -deps tvdenoisefig1.eps
%
% Add noise.
%
noiselevel=50;
img2=img+noiselevel*randn(rows,cols);
%
% Round back to 8 bit gray levels from 0 to 255.
%
img2=round(img2);
img2(img2<0)=0;
img2(img2>255)=255;
%
% Plot the noisy image.
%

figure(2);
clf
imagesc(img2);
colormap gray;
bookfonts
axis off

disp('Displaying the noisy image (fig. 2)')
print -deps tvdenoisefig2.eps
%
% Setup the denoising problem.  The data is the noisy image.  The G
% matrix is simply the identity (in sparse form.)  The L matrix
% does finite difference approximations to the derivatives in the
% horizontal and vertical directions.  
%
dn=reshape(img2,rows*cols,1);
G=speye(rows*cols);
%
% The L matrix is large and sparse but somewhat complicated.  One
% way to build it is using the sparse() function and supplying the
% i,j indices and entry values as vectors.
%
% There are (row-1)*(cols-1)*4 nonzero entries. So we allocate this much space.
%
numentries=(rows-1)*(cols-1)*4;
iindices=zeros(numentries,1);
jindices=zeros(numentries,1);
entries=zeros(numentries,1);
%
% Start with the first row of the matrix and work our way up.
%
r=1;
%
% Keep track of the number of entries.
%
k=1;
%
% Setup to use sub2ind
%
sz=size(img);
%
% Loop over the pixels.
%
for j=1:cols-1
  for i=1:rows-1
    %
    % m(i,j)
    %
    iindices(k)=r;
    jindices(k)=sub2ind(sz,i,j);
    entries(k)=1;
    k=k+1;
    %
    % -m(i,j+1)
    %
    iindices(k)=r;
    jindices(k)=sub2ind(sz,i,j+1);
    entries(k)=-1;
    k=k+1;
    %
    % Next row, m(i,j)-m(i+1,j)
    %
    r=r+1;
    %
    % m(i,j)
    %
    iindices(k)=r;
    jindices(k)=sub2ind(sz,i,j);
    entries(k)=1;
    k=k+1;
    %
    % -m(i+1,j)
    %
    iindices(k)=r;
    jindices(k)=sub2ind(sz,i+1,j);
    entries(k)=-1;
    k=k+1;
    %
    % Next row.
    %
    r=r+1;
  end
end
%
% Add an entry to access the bottom right corner.
%
iindices(k)=r;
jindices(k)=rows*cols;
entries(k)=0.0;
%
% Now, build the L matrix for TV regularization
%
L=sparse(iindices,jindices,entries);

%alpha loop to examine the L-curve for this problem, if desired
% alphas = logspace(0,3,8);
% for i=1:length(alphas)
%     [mtmp,bestobj,iter]=admml1reg(G,dn,L,alphas(i));
%     misfit(i)=norm(G*mtmp-dn);
%     onenorm(i)=norm(L*mtmp,1);
% end
% 

% figure(3)
% clf
% plot(misfit,onenorm,'*');
% ylabel('|| Lm ||_1')
% xlabel('|| Gm - d ||');
% bookfonts
    
%
% Try alpha=3.0;  Too Small
%
[mreg03,bestobj03,iter03]=admml1reg(G,dn,L,2.5);
img03=round(reshape(mreg03,512,512));
img03(img03<0)=0;
img03(img03>255)=255;

disp('Displaying the denoised image for alpha=2.5 (fig. 10)')

figure(4);
clf
imagesc(img03);
colormap(gray);
bookfonts
axis off;
print -deps tvdenoisefig10.eps
%
% Try alpha=50. Nearly optimal.
%
[mreg50,bestobj50,iter50]=admml1reg(G,dn,L,50);
img50=round(reshape(mreg50,512,512));
img50(img50<0)=0;
img50(img50>255)=255;

disp('Displaying the denoised image for alpha=50 (fig. 11)')

figure(5);
clf
imagesc(img50);
bookfonts
colormap(gray);
axis off
print -deps tvdenoisefig11.eps

%
% Try alpha=300. Too big.
%
[mreg300,bestobj300,iter300]=admml1reg(G,dn,L,300);
img300=round(reshape(mreg300,512,512));
img300(img300<0)=0;
img300(img300>255)=255;

disp('Displaying the denoised image for alpha=300 (fig. 12)')

figure(6);
clf
imagesc(img300);
colormap(gray);
axis off
print -deps tvdenoisefig12.eps
