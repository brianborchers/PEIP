% Example 6.3 part a
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% This contains the system setup and data creation for a large-scale 
% tomography problem

% make sure we have a clean environment
clear
rand('seed',0);
randn('seed',0);

% block edges for this problem are defined at integer boundaries 
% between 0 and msize
msize=30;
noise=0.05;
%block center coordinates

disp('Generating blob model')

%initialize the model
mtrue=zeros(msize,msize,msize);

%modify the background to get a random blob model
var=3;
mblob=100;
for l=1:mblob
  ir1=randn*3;
  ir2=randn*3;
  ir3=randn*3;
  a=randn;
  for i=1:msize
    for j=1:msize
      for k=1:msize
        mtrue(i,j,k)=mtrue(i,j,k)+...
          a*exp(...
          (-(i-ir1-msize/2)^2-(0.5*(j-ir2-msize/2))^2-(k-ir3-msize/2)^2)/var);
      end
    end
  end
end

%make the model have zero mean
mtrue=mtrue-mean(mean(mean(mtrue)));

%reshape it into a column vector
mtrue_rs=reshape(mtrue,msize^3,1);


%calculate rays
%start and end points in 3-d are stored in columns of arrays x1 and x2

%rays in each rayset (there are three (x,y, and z-spanning) raysets

disp('Generating random raypaths');
N=10000;
%random rays beginning and ending at block centers at edge of model
%full z spanning
x11z=floor(msize*rand(1,N))+rand(1,N);
x21z=floor(msize*rand(1,N))+rand(1,N);

x12z=floor(msize*rand(1,N))+rand(1,N);
x22z=floor(msize*rand(1,N))+rand(1,N);

x13z=zeros(1,N);
x23z=msize*ones(1,N);

%full y spanning
x11y=floor(msize*rand(1,N))+rand(1,N);
x21y=floor(msize*rand(1,N))+rand(1,N);

x12y=zeros(1,N);
x22y=msize*ones(1,N);

x13y=floor(msize*rand(1,N))+rand(1,N);
x23y=floor(msize*rand(1,N))+rand(1,N);

%full x spanning
x11x=zeros(1,N);
x21x=msize*ones(1,N);

x12x=floor(msize*rand(1,N))+rand(1,N);
x22x=floor(msize*rand(1,N))+rand(1,N);

x13x=floor(msize*rand(1,N))+rand(1,N);
x23x=floor(msize*rand(1,N))+rand(1,N);

% put everything into the main raypath vectors
x1=[[x11z;x12z;x13z],[x11y;x12y;x13y],[x11x;x12x;x13x]];
x2=[[x21z;x22z;x23z],[x21y;x22y;x23y],[x21x;x22x;x23x]];

[~,nrays]=size(x1);

%matrix of ray vectors between x1 and x2
v=x2-x1;

disp('Populating the G matrix and travel time vector; this may takes a while')
%initilize G matrix and travel time vector
GT=sparse(msize^3,nrays);
tt=zeros(nrays,1);

%populate G matrix
for ray=1:nrays
  % first determine the noiseless travel time
  % this is done by determining the length in each cell each ray traverses
  % then multiplying by the block slowness
  
  % raytrace by finding exact solutions for integer block boundaries 
  % in all three dimensions
  % alpha holds the coefficients where x1+alpha*v contains at least one integer
  for dimension=1:3
    for d=0:msize
      alpha(dimension,d+1)=(d-x1(dimension,ray))/v(dimension,ray);
    end
  end

  % we only care about the unique alphas, and will want them in order
  alphavec=unique(sort(reshape(alpha,3*(msize+1),1)));

  % select the appropriate range of alpha values to ray trace
  % this removes values for cells that weren't intersected by the ray
  % they are introduced because we search all boundaries when only some occur
  % within the scanned region
  ind=find(alphavec >=0 & alphavec <=1);
  alphavec=alphavec(ind);

  % segment lengths for insertion into G matrix
  lengths=alphavec*norm(v(:,ray));
  segs=diff(lengths);

  % get the block each segment is in by getting the midpoint 
  % midpoints for each ray segment and taking the ceiling of it
  dalphavec=diff(alphavec);
  for i=1:length(alphavec)-1
    midpoint(:,i)=x1(:,ray)+(alphavec(i+1)-dalphavec(i)/2)*v(:,ray);
  end
  blocks=ceil(midpoint);

  % now generate the 3d model version of the model for this ray
  G3r=zeros(msize,msize,msize);
  for i=1:length(segs)
    G3r(blocks(1,i),blocks(2,i),blocks(3,i))=segs(i);
  end

  % add an appropriate column to the GT matrix.
  Grow=reshape(G3r,msize^3,1);
  GT(:,ray)=Grow';

  % compute the true travel time
  tt(ray)=sum(sum(sum(mtrue.*G3r)));
end

% Get G from GT.  
G=GT';

%add data noise
ttn=tt+noise*randn(size(tt));

% We'll use "d" for the right hand side of the least squares problem.
d=ttn;

disp('creating the regularization matrix');
% Make an L matrix with Laplacian smooth and a bit of 0th order regularization.
Lrows=(msize-2)^3+msize^3;

% Make the LT matrix, then transpose to get L.
LT=sparse(msize^3,Lrows);

% First add the Laplacian smoothing.
zerorow=zeros(msize,msize,msize);
currentrow=1;
for i=2:msize-1
  for j=2:msize-1
    for k=2:msize-1
      row=zerorow;
      row(i-1,j,k)=1;
      row(i+1,j,k)=1;
      row(i,j-1,k)=1;
      row(i,j+1,k)=1;
      row(i,j,k-1)=1;
      row(i,j,k+1)=1;
      row(i,j,k)=-6;
      LT(:,currentrow)=reshape(row,msize^3,1);
      currentrow=currentrow+1;
    end
  end
end

% now, the 0th order regularization.
for i=1:msize^3
  LT(i,currentrow)=0.01;
  currentrow=currentrow+1;
end

% Now, transpose to get L.
L=LT';

% save what will be needed for the next example script (ex_6_4b)
save ex_6_4_30000.mat G L tt mtrue d
