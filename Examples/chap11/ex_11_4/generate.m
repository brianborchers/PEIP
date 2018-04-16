% Example 11.4 proposal generator
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% y=generate(x)
%
%
% For this problem, we'll use a multivariate normal generator, with 
% standard deviations specified by the vector step.

% Note that logproposal.m and generate.m
% are closely tied to each other.
%
function y=generate(x)
global step
y=x+step.*randn(4,1);
