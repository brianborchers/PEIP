% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% X=idcost(Y)
%
%   Takes the column-by-column inverse discrete cosine transform of Y
%
function X = idcost(Y)

[n, m] = size(Y);
X = zeros(size(Y));

% Precompute weights
ww = sqrt(2*n) * exp(1i*(0:n - 1)*pi/(2 * n)).';

%loop over columns
for i = 1:m
    if (rem(n, 2) == 1)
        % odd case
        if i == 1, ww(1) = ww(1) * sqrt(2); end
        yy = zeros(2*n, 1);
        yy(1:n, :) = ww .* Y(:, i);
        yy(n+2:2*n, :) = - 1i * ww(2:n) .* Y(n:-1:2, i);
        
        y = ifft(yy);
    else
        % even case
        if i == 1, ww(1) = ww(1) / sqrt(2); end
        yy = ww .* Y(:, i);
        yt = ifft(yy);
        y = zeros(n, 1);
        y(1:2:n) = yt(1:n/2);
        y(2:2:n) = yt(n:-1:n/2+1);
    end
    
    % Extract inverse DCT getting rid of the imaginary roundoff error
    X(:, i) = real(y(1:n));
end
