% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% y=admml1regmult(x,mode)
%
% Computes y=A*x or y=A'*x for the Alternating Direction Method of
% Multipliers algorithm.  If mode=='transp', then A'*x is used.  If
% mode=='notransp' then A*x is used.
%
% In terms of G, L, R, and alpha,
%
% A=[G; sqrt(rho/2)*L]
%
% So,
%
% A*x=[G*x; sqrt(rho/2)*L*x]
%
% A'*x=[G'*x(1:n) + sqrt(rho/2)*(L'*x(n+1:end))]
%
function y = admml1regmult(x, mode)
%
% Global variables hold G, L, W, and ALPHA.
%
global ADMML1REGG;
global ADMML1REGL;
global ADMML1REGRHO;
global ADMML1REGGT;
global ADMML1REGLT;
%
%
% Check to see which mode we're in and do the right thing.
%
switch mode
    %
    % First case, no transpose.
    %
    case 'notransp'
        y = [ADMML1REGG * x; sqrt(ADMML1REGRHO/2) * (ADMML1REGL * x)];
        return
        %
        % Second case, we've got a transpose. We've precomputed the
        % transposes of G and L and padded with 0's, so that we can just do
        % [G' 0]*x + sqrt(rho/2)*[0 L']*x
        %
    case 'transp'
        y = ADMML1REGGT * x + sqrt(ADMML1REGRHO/2) * (ADMML1REGLT * x);
        return
end
