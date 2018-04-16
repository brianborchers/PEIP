% Example 6.1, 
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber

% make sure we have a clean environment
clear
rand('state',0);
randn('state',0);

% Number of model parameters in the problem
m = 256;

% Initialize the design matrix
% Parameter indices increase going down each column
% G(:,1) is in the NW corner, G(:,16) is in the SW corner
% G(:,241) in the NE corner, and G(:,256) is in the SE corner
G = zeros(94,m);

% the row and columns scans have 1 where a ray runs through a block
% Design matrix entries for the column scan
for i=1:16
  for j=(i-1)*16+1:i*16
    G(i,j) = 1.;
  end
end

% Design matrix for the row scan
for i=1:16
  for j=i:16:240+i
    G(i+16,j) = 1.;
  end
end

% The diagonal scans, because these are going through on the diagonal the 
% entries are sqrt(2)

% G matrix for the SW to NE diagonal scan, upper part
for i=1:16
    for j=0:i-1
        G(i+32,i+j*15) = sqrt(2.);
    end
end

% G matrix for the SW to NE diagonal scan, lower part
for i=1:15
    for j=0:15-i
        G(i+48,(i+1)*16+j*15) = sqrt(2.);
    end
end

% G matrix for the NW to SE diagonal scan, lower part
for i=1:16
    for j=0:i-1
        G(i+63,17-i+17*j) = sqrt(2.);
    end
end

% G matrix for the NW to SE diagonal scan, upper part
for i=1:15
    for j=0:15-i
        G(i+79,(i*16)+1+17*j) = sqrt(2.);
    end
end

% Setup the true model.
mtruem=zeros(16,16);
mtruem(9,9)=1;
mtruem(9,10)=1;
mtruem(9,11)=1;
mtruem(10,9)=1;
mtruem(10,11)=1;
mtruem(11,9)=1;
mtruem(11,10)=1;
mtruem(11,11)=1;
mtruem(2,3)=1;
mtruem(2,4)=1;
mtruem(3,3)=1;
mtruem(3,4)=1;

% reshape the true model to be a vector
mtruev=reshape(mtruem,256,1);

% Compute the data.
dtrue=G*mtruev;

% Add the noise.
d=dtrue+0.01*randn(size(dtrue));

% Plot the true model.
figure(1)
clf
bookfonts
imagesc(mtruem,[-0.2 1.2]);
colormap(gray);
H=colorbar;
set(H,'FontSize',18);

disp('Displaying true model (fig. 1)')
print -deps mtomo.eps

% Next, compute the generalized inverse solution, using 87 singular values.

% t is starting cputime for finding computation time 
t=cputime;

% Get svd
[U,S,V]=svd(G);

% Truncate to 87 singular values
UP=U(:,1:87);
VP=V(:,1:87);
SP=S(1:87,1:87);

% Get solution using 87 singular values
mg=VP*inv(SP)*UP'*d;

% Show computation time for gen inv solution
disp(['cputime for gen inv solution was ', num2str(cputime-t)]);

% And plot the solution.
figure(2)
clf
imagesc(reshape(mg,16,16),[-0.2 1.2]);
bookfonts
colormap(gray);
H=colorbar;
set(H,'FontSize',18);
disp('Displaying truncated svd solution with 87 singular values (fig. 2)')
print -deps mg.eps

% Next,compute the Kaczmarz solution, 94 iterations.

% t is starting cputime for finding computation time 
t=cputime;

% Get Kaczmarz solution, the 0.0 tolerance forces 200 iterations to happen
mkac=kac(G,d,0.0,200);

% Show computation time for Kaczmarz solution
disp(['cputime for Kaczmarz''s algorithm was ', num2str(cputime-t)]);

% And plot the solution.
figure(3)
clf
imagesc(reshape(mkac,16,16),[-0.2 1.2]);
bookfonts
colormap(gray);
H=colorbar;
set(H,'FontSize',18);
disp('Displaying Kaczmarz algorithm solution (fig. 3)')
print -deps mkac.eps

% Display relative errors for different solution methods
disp(' ')
disp(['kac relative error  ', num2str(norm(mkac-mtruev)/norm(mtruev))])
disp(['SVD relative error  ', num2str(norm(mg-mtruev)/norm(mtruev))])
