% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% plotconst(x,l,r)
%
% Plots a model in piecewise constant form over n subintervals,
% where n is the length of x.
%
% Input Parameters:
%   x - model to be plotted
%   l - left endpoint of plot
%   r - right endpoint of plot
%
% Output Parameters: None.  This function returns no value.
function H = plotconst(x, l, r)

% Find length of model
n = length(x);
% Find size of each interval
delta = (r - l) / n;
% Dummy values at beginning of vectors to allow concatination
myx = [0];
myy = [0];
% Iteratively fill vector of x and y values for steps for plot
% by concatenating onto dummy vectors
for i = 1:n
    myx = [myx,((i - 1)*delta+l:(delta / 20):i*delta+l)];
    myy = [myy,(ones(1, 21)*x(i))];
end

% Find length of resulting vector of x values
l2 = length(myx);
% Truncate vectors to remove dummy values used in concatination
myx = myx(2:l2);
myy = myy(2:l2);

% Plot piecewise constant graph
H = plot(myx, myy, 'k-');
