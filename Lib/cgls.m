% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% [X,rho,eta]=cgls(G,d,niter)
%
% Performs niter iterations of the CGLS algorithm on the least
% squares problem
%
%   min norm(G*m-d)
%
% The iterates 1, 2, ..., niter are returned in the columns of the
% matrix X.  For each iterate we also compute rho(i)=norm(G*m-d)
% and eta(i)=norm(m).
function [X, rho, eta] = cgls(G, d, niter)
%

% Figure out problem size.
[nrows, ncols] = size(G);
if (length(d) ~= nrows)
    error('G and d do not match in size.');
end

% Setup space for the results.
X = zeros(ncols, niter);
rho = zeros(niter, 1);
eta = zeros(niter, 1);

% Setup for the first iteration.
m = zeros(ncols, 1);
p = zeros(ncols, 1);
beta = 0;
s = d;
r = G' * s;

% Main loop- perform CGLS iterations.
for k = 0:niter - 1
    % We'll precompute r'*r since it's used in several places.
    rtr = r' * r;
    
    %  Update beta.
    if (k > 0)
        beta = rtr / (prevr' * prevr);
    end
    
    %  Update p
    p = r + beta * p;
    
    % Compute the new alpha.  To avoid doing the matrix vector
    % multiplication repeatedly, we store G*p in Gp
    Gp = G * p;
    alpha = rtr / (Gp' * Gp);
    
    % Update m.
    m = m + alpha * p;
    
    % Update s.
    s = s - alpha * Gp;
    
    % Save r for the next iteration, and then update it.
    prevr = r;
    r = G' * s;
    
    % Store the new iterate.
    X(:, k+1) = m;
    rho(k+1) = norm(s);
    eta(k+1) = norm(m);
end
