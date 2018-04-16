% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% y=irlsl1regmult(x,mode)
%
% Computes y=A*x or y=A'*x for the iteratively reweighted least
% squares algorithm.  If mode=='transp', then A'*x is used.  If
% mode=='notransp' then A*x is used.
%
% In terms of G, L, R, and alpha,
%
% A=[G; sqrt(alpha)*sqrt(R)*L]
%
% So,
%
% A*x=[G*x; sqrt(alpha)*sqrt(R)*L*x]
%
% A'*x=[G'*x+\sqrt(alpha)*L'*(W*x)]
%
% Note that W is a diagonal weighting matrix stored as a simple
% vector,  so we can use .* rather than full matrix
% multiplication for those multiplications involving W. Also note
% that we've folded the factor of sqrt(alpha/2) into W to make this
% run slightly faster.
%
%
function y = irlsl1regmult(x, mode)
%
% Global variables hold G, L, W, and ALPHA.
%
global IRLSL1REGG;
global IRLSL1REGL;
global IRLSL1REGW;
%
%
% Check to see which mode we're in and do the right thing.
%
switch mode
    %
    % First case, no transpose.
    %
    case 'notransp'
        y = [IRLSL1REGG * x; IRLSL1REGW .* (IRLSL1REGL * x)];
        return
        %
        % Second case, we've got a transpose.  We need to break x up into
        % two parts.
        %
    case 'transp'
        n = size(IRLSL1REGG, 1);
        y = IRLSL1REGG' * x(1:n) + IRLSL1REGL' * (IRLSL1REGW .* x((n + 1):end));
        return
end
