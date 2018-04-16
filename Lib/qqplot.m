% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% qqplot(data)
%
% Produces a Q-Q plot of a data vector compared with a standard normal distribution
function qqplot(data)
m = mean(data);
s = std(data);
n = length(data);
x = sort(data);
q = zeros(n, 1);
xqm = zeros(n, 1);
for i = 1:length(x)
    q(i) = (i - 0.5) / n;
    %compare against a normal distribution
    xqm(i) = phiinv(q(i));
end
plot(xqm, x, 'k+');
hold on
%
% Plot the theoretical line for the given mean and standard deviation
%
plot([xqm(1), xqm(n)], [m + phiinv(0.5/n) * s, m + phiinv((n - 0.5)/n) * s], 'k--');

