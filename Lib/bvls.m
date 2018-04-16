% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% x=bvls(A,b,l,u,maxiter)
%
% Solves a bounded variables least squares problem
%
%   min || Ax-b ||
%          l <= x <= u
%
% The algorithm is based on the bounded variables least squares algorithm
% of Stark and Parker.  See
%
%  P.B. Stark and R. L. Parker, "Bounded-Variable Least-Squares: An Algorithm
%  and Applications", Computational Statistics 10:129-141, 1995.
%
% INPUT
%   A       - a system matrix
%   b       - the right hand side
%   l       - the lower bounds on x
%   u       - the upper bounds on x
%   maxiter - the maximum number of iterations to perform
%
% OUTPUT
%   x     - the solution found
%   state - 0 if the corresponding x is between its bounds, 1 if at its lower
%           bound and 2 if at its upper bound
%
function [x, state] = bvls(A, b, l, u, maxiter)

% Get the size of A for future reference.
[~, n] = size(A);

% If maxiter isn't specified, use a default value.
if (nargin < 5)
    maxiter = 10 * n;
end

% a tolerance setting
myeps = 1.0e-10;

% Initialize a bunch of variables.  This speeds things up by avoiding
% reallocation of storage each time an entry is added to the end of a vector.
state = zeros(n, 1);
x = zeros(n, 1);
criti = 0;

% vectors of indices that we have to grow
oopslist = [];
atbound = [];
between = [];

% setup an initial solution with all vars locked at a lower or upper bound.
% set any free vars to 0.
for i = 1:n
    % determine the starting value for x(i)
    if ((u(i) >= + Inf) & (l(i) <= - Inf))
        % we are a free variable
        x(i) = 0;
        state(i) = 0;
        between = [between, i];
    elseif (abs(l(i)) < abs(u(i)))
        % if the lower bound is closer 0 than the upper bound use it
        x(i) = l(i);
        state(i) = 1;
        atbound = [atbound, i];
    else
        % otherwise use the upper bound is no further than the lower bound so use it
        x(i) = u(i);
        state(i) = 2;
        atbound = [atbound, i];
    end
end

% The main loop.  Stop after maxiter iterations if nothing else.
iter = 0;
while (iter < maxiter)
    iter = iter + 1;
    
    % optimality test.
    grad = A' * (A * x - b);
    
    % We ignore any variable that has already been tried and failed.
    grad(oopslist) = zeros(size(oopslist));
    
    % Check for optimality.
    done = 1;
    for i = 1:n
        % if this element of x is free and has noticeable nonzero gradient
        if ((abs(grad(i)) > ((1 + norm(b)) * myeps)) & (state(i) == 0))
            done = 0;
            break;
        end
        % if we are at our lower bound and have negative gradient
        if ((grad(i) < 0) & (state(i) == 1))
            done = 0;
            break;
        end
        % if we are at our upper bound and have negative gradient
        if ((grad(i) > 0) & (state(i) == 2))
            done = 0;
            break;
        end
    end
    
    % if no indices were invalid we have converged
    if (done == 1)
        return;
    end
    
    % Not optimal, so we need to free up some vars at bounds.
    newi = 0;
    newg = 0.0;
    % Look for the locked variable with the biggest gradient.
    for i = atbound
        if ((grad(i) > 0) & (state(i) == 2))
            if (abs(grad(i)) > newg)
                newi = i;
                newg = abs(grad(i));
            end
        end
        if ((grad(i) < 0) & (state(i) == 1))
            if (abs(grad(i)) > newg)
                newi = i;
                newg = abs(grad(i));
            end
        end
    end
    
    % Free the locked variable with the biggest gradient if there is one.
    if (newi ~= 0)
        atbound = remove(atbound, newi);
        state(newi) = 0;
        between = [between, newi];
    end
    
    % Make sure the projected problem is nontrivial.
    if (length(between) == 0)
        disp('Empty projected problem');
        continue;
    end
    
    % Construct the new projected problem.
    Aproj = A(:, between);
    An = A(:, atbound);
    if (length(atbound) > 0)
        bproj = b - An * x(atbound);
    else
        bproj = b;
    end
    
    % Solve the projected problem.
    z = Aproj \ bproj;
    
    % populate the new x
    xnew = zeros(size(x));
    xnew(atbound) = x(atbound);
    xnew(between) = z;
    
    % if there was an unlocked variable and after solving the projected
    % problem it went beyond the bound it was released from it moved an illegal
    % direction so needs to be blocked until we find a variable that will move
    % legally
    if ((newi ~= 0) & ...
            (((xnew(newi) <= l(newi)) & (x(newi) == l(newi))) | ...
            ((xnew(newi) >= u(newi)) & (x(newi) == u(newi)))))
        
        oopslist = [oopslist; newi];
        if ((xnew(newi) <= l(newi)) & (state(newi) == 1))
            state(newi) = 1;
            x(newi) = l(newi);
        end
        if ((xnew(newi) >= u(newi)) & (state(newi) == 2))
            state(newi) = 2;
            x(newi) = u(newi);
        end
        atbound = [atbound, newi];
        between = remove(between, newi);
        
        % go back to the top of the while loop
        continue;
    end
    
    % We've got a good variable freed up.  Reset the oopslist.
    oopslist = [];
    
    % Move as far as possible towards the optimal solution to the projected
    % problem.
    
    % alpha is the farthest fraction of the distance moved that keeps x in bounds
    alpha = 1;
    for i = between
        % if we were above our upper bound
        if (xnew(i) > u(i))
            newalpha = min([alpha, (u(i) - x(i)) / (xnew(i) - x(i))]);
            if (newalpha < alpha)
                criti = i;
                crits = 2;
                alpha = newalpha;
            end
        end
        if (xnew(i) < l(i))
            newalpha = min([alpha, (l(i) - x(i)) / (xnew(i) - x(i))]);
            if (newalpha < alpha)
                criti = i;
                crits = 1;
                alpha = newalpha;
            end
        end
    end
    
    % Take the step.
    x = x + alpha * (xnew - x);
    
    % Update the state of limiting index
    if (alpha < 1)
        between = remove(between, criti);
        atbound = [atbound, criti];
        state(criti) = crits;
    end
    
    % Update the state (and value if necessary) of all indices
    for i = 1:n
        if (x(i) >= u(i))
            x(i) = u(i);
            state(i) = 2;
            
            % ensure that i is only in atbound
            if (~isempty(find(between == i)))
                between = remove(between, i);
            end
            if (isempty(find(atbound == i)))
                atbound = [atbound, i];
            end
        end
        if (x(i) <= l(i))
            x(i) = l(i);
            state(i) = 1;
            
            % ensure that i is only in atbound
            if (~isempty(find(between == i)))
                between = remove(between, i);
            end
            if (isempty(find(atbound == i)))
                atbound = [atbound, i];
            end
        end
    end
    
    % continue iterating
end
disp('BVLS Exceeded max iters')
