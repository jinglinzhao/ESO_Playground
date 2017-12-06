% This function is not recommended because the output seems to be noisier
% than the FFT_noise function. 
% Tested with Gaussian and cosine funciton. 
% https://au.mathworks.com/help/matlab/math/fourier-transforms.html
% https://au.mathworks.com/help/matlab/ref/fft.html

function [f,P] = FFT_g(X, delta_t)

L   = length(X);        % Signal length
Fs = 1 / delta_t;       % Sampling frequency

% To use the fft function to convert the signal to the frequency domain, 
% first identify a new input length that is the next power of 2 from the 
% original signal length. This will pad the signal X with trailing zeros 
% in order to improve the performance of fft.
n = 2^nextpow2(L);  % can be tuned 
Y = fft(X,n);           % n-point discrete Fourier transform of X
f = Fs*(0:(n/2))/n;     % Frequency 
P = abs(Y/n);           % Power
P = P(1:n/2+1);

end