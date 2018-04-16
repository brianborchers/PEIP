% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%Generate the m by n system matrix for the Shaw problem
%
% G = shaw_G(n,m)
%
%
% m and n must be even.

% Reference: C. B. Shaw, Jr., "Improvements of the resolution of
% an instrument by numerical solution of an integral equation",
% J. Math. Anal. Appl. 37m 83-112, 1972.
function G = shaw_G(m, n)

if (rem(m, 2) ~= 0), error('m must be even'), end
if (rem(n, 2) ~= 0), error('n must be even'), end

% Initialization.
delta = pi / n;
G = zeros(m, n);
vm = ((1:m) - 1 / 2) * pi / m - pi / 2;
vn = ((1:n) - 1 / 2) * pi / n - pi / 2;

for i = 1:m
    for j = 1:n
        s = vm(i);
        th = vn(j);
        x = pi * (sin(s) + sin(th));
        if x == 0
            G(i, j) = (cos(s) + cos(th))^2;
        else
            G(i, j) = (cos(s) + cos(th))^2 * (sin(x) / x)^2;
        end
    end
end
G = G * delta;
