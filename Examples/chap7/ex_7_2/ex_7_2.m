% Example 7.3
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber

% make sure we have a clean environment
clear
rand('state',0);
randn('state',0);

% Load the image.
img=double(imread('image.png'));

% Display the image.
figure(1)
clf
imagesc(img);
colormap gray;
axis off
bookfonts

disp('Displaying the original image (fig. 1)')
print -depsc c7foriginalimg.eps

% Perform discrete cosine transform threshholding.
% Here, we zero out the smallest 40% of the DCT coefficients
% and show that it has only a very small effect on the image.
img2=dctthresh(img,0.4);

% Display the thresholded image
figure(2)
clf
imagesc(img2)
bookfonts
colormap gray;
axis off

disp('Displaying the sparser image (fig. 2)')
print -depsc c7fthreshimg.eps
