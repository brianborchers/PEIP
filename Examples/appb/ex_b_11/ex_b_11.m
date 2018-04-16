% Script file for PEIP Example B.11
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% Requires MATLAB Statistics Toolbox
%Reinitialize everything.
%
clear
%rand('seed',0);
%randn('seed',0);
%
% Generate the data set.
%
%Student's T distribution data, 1000 points
%Requires Statistics Toolbox
%data=trnd(5,1000,1);
load data
%
% Plot the histogram.
%
figure(1);
clf
bookfonts;
disp('Displaying Histogram of Sample Data (fig 1)')
hist(data,(-8:8))
h = findobj(gca,'Type','patch');
set(h,'FaceColor','k')
ylabel('N')
xlabel('x')
bookfonts

print -deps2 abfhist.eps

figure(2);
clf
bookfonts;
qqplot(data);
axis tight
disp('Displaying Q-Q Plot of Sample Data vs. Standard Normal (fig 2)')
xlabel('Standard Normal Quantiles')
ylabel('Quantiles of Input Sample')
bookfonts
print -deps2 abfQQ.eps



