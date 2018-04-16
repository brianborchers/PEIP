% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% w=bartl(m)
%
%returns the Bartlett (triangle) window of length m
%
function w = bartl(m)
w = 2 * (0:(m - 1) / 2) / (m - 1);
if rem(m, 2) == 1
    w = [w, w((m - 1)/2:-1:1)]';
else
    w = [w, w(m/2:-1:1)]';
end
