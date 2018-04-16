% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% plotconstc(x,l,r,c)
%
% Plots a model in piecewise constant form over n subintervals,
% where n is the length of x.  Line type and color for plotting
% are an input parameter.
%
% This also sets the LineWidth to 1 and the FontSize to 18
%
% Input Parameters:
%   x - model to be plotted
%   l - left endpoint of plot
%   r - right endpoint of plot
%   c - line type and color parameter
%
% Output Parameters:
%   H - handle for the plot.
%
function H = plotconstc(x, l, r, c)

% Find length of model

n = length(x);

% Find size of each interval
delta = (r - l) / n;

% Dummy values at beginning of vectors to allow concatination
myx = [0];
myy = [0];

% Iteratively fill vector of x and y values for steps for plot
% by concatenating the existing vectors for each value add 21 points at the
% same height
for i = 1:n
    myx = [myx,((i - 1)*delta+l:(delta / 20):i*delta+l)];
    myy = [myy,(ones(1, 21)*x(i))];
end

% Remove the 0 elements from the beginning added to allow concatenation
myx = myx(2:end);
myy = myy(2:end);

% Plot piecewise constant graph with input line type and color
plot(myx, myy, c);
H = gca;
set(H, 'FontSize', 18);
set(H, 'LineWidth', 1.0);
