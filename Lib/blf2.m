% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% [xmin,xmax]=blf2(A,b,c,delta,l,u)
%
% Solves the problem
%
%   min/max c'*x
%    || Ax-b || <= delta
%    l <= x <= u
%
% Based on an algorithm described in
%
%
%  P.B. Stark and R. L. Parker, "Bounded-Variable Least-Squares: An Algorithm
%  and Applications", Computational Statistics 10:129-141, 1995.
%
function [xmin, xmax] = blf2(A, b, c, delta, l, u)
x0 = bvls(A, b, l, u);
gamma0 = c' * x0;
alpha = 1.0e-3;

% First, make sure that || A x0 - b || <= delta.  If not, then the problem
% is just infeasible.
if (norm(A*x0-b) > delta)
    xmin = [];
    xmax = [];
    return;
end


% now find xmax
% Find a gamma so that || A*x(gamma)-b || > delta
step = 0.1 * abs(gamma0) + 0.01;
gamma = gamma0 + step;
G = [A; alpha * c'];
d = [b; alpha * gamma];
xgamma = bvls(G, d, l, u);
xgammanorm = norm(A*xgamma-b);

% left is a step added to gamma known to have xgammanorm < delta
left = 0;
while ((xgammanorm < delta) & (step < 1.0e10))
    left = step;
    step = step * 2;
    gamma = gamma0 + step;
    G = [A; alpha * c'];
    d = [b; alpha * gamma];
    xgamma = bvls(G, d, l, u);
    xgammanorm = norm(A*xgamma-b);
end
right = step;

%  Now, do a binary search to narrow the range.
while ((right - left) / (0.0001 + right) > 0.01)
    mid = (left + right) / 2;
    gamma = gamma0 + mid;
    G = [A; alpha * c'];
    d = [b; alpha * gamma];
    xgamma = bvls(G, d, l, u);
    xgammanorm = norm(A*xgamma-b);
    if (xgammanorm > delta)
        right = mid;
    else
        left = mid;
    end
end

% Save this solution as xmax
gamma = gamma0 + left;
G = [A; alpha * c'];
d = [b; alpha * gamma];
xmax = bvls(G, d, l, u);
cmax = c' * xmax;


% Now find xmin
% Find a gamma so that || A*x(gamma)-b || > delta
%
step = 0.1 * abs(gamma0) + 0.01;
gamma = gamma0 - step;
G = [A; alpha * c'];
d = [b; alpha * gamma];
xgamma = bvls(G, d, l, u);
xgammanorm = norm(A*xgamma-b);

% left is a step added to gamma known to have xgammanorm < delta
left = 0;
while ((xgammanorm < delta) & (step < 1.0e10))
    left = step;
    step = step * 2;
    gamma = gamma0 - step;
    G = [A; alpha * c'];
    d = [b; alpha * gamma];
    xgamma = bvls(G, d, l, u);
    xgammanorm = norm(A*xgamma-b);
end
right = step;

%  Now, do a binary search to narrow the range.
while ((right - left) / (0.0001 + right) > 0.01)
    mid = (left + right) / 2;
    gamma = gamma0 - mid;
    G = [A; alpha * c'];
    d = [b; alpha * gamma];
    xgamma = bvls(G, d, l, u);
    xgammanorm = norm(A*xgamma-b);
    if (xgammanorm > delta)
        right = mid;
    else
        left = mid;
    end
end

% Save this solution as xmin
gamma = gamma0 - left;
G = [A; alpha * c'];
d = [b; alpha * gamma];
xmin = bvls(G, d, l, u);
cmin = c' * xmin;


% ensure that xmin corresponded with the lower dot average
if (cmin > cmax)
    temp = xmax;
    xmax = xmin;
    xmin = temp;
end
