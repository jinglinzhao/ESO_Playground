% Hermite function normalized to 1. 
% Width is not included in the parameter

function y = hermite_nor(n, x)
    y = (sqrt(pi) * 2^n * factorial(n)).^(-0.5) .* hermiteH(n, x) .* exp(-x.^2/2);
end