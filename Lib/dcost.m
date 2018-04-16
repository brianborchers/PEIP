% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%  Y=dcost(X)
%
%   computes the column-by-column discrete cosine transform of X.
%
function Y = dcost(X)

[n, m] = size(X);
Y = zeros(size(X));

% Precompute weights
ww = (exp(-1i*(0:n - 1)*pi/(2 * n)) / sqrt(2*n)).';
ww(1) = ww(1) / sqrt(2);

%loop over columns
for i = 1:m
    if (rem(n, 2) == 1)
        % odd case
        y = zeros(2*n, 1);
        y(1:n) = X(:, i);
        y(n+1:2*n) = X(n:-1:1, i);
        ff = fft(y);
        ff = ff(1:n);
    else
        % even case
        y = [X(1:2:n, i); X(n:-2:2, i)];
        ff = fft(y);
        if i == 1, ww = 2 * ww; end
    end
    
    Y(:, i) = real(ww.*ff);
end

