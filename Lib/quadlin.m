% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% Function [x,lambda]=quadlin(Q,A,b) solves the problem:
%   min (1/2) x'*Q*x
%              Ax=b
% using the Lagrange multiplier technique, where Q is assumed to be
% symmetric and positive definite and the rows of A are linearly
% independent.
%
% Input Parameters:
%   Q - positive definite symmetric matrix
%   A - matrix with linearly independent rows
%   b - data vector
%
% Output Parameters:
%   x - vector of solution values
%   lambda - Lagrange multiplier
%
% The formulas for the solution of this problem are:
%
function [x, lambda] = quadlin(Q, A, b)

% First, find lambda.
lambda = (A * inv(Q) * A') \ b;

% Now, x.
x = inv(Q) * A' * lambda;
