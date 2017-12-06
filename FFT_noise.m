function [fshift,power] = FFT_noise(xnoise, delta_t)

Fs      = 1/delta_t;
n       = length(xnoise);
fshift  = (-n/2:n/2-1)*(Fs/n);

% Power is the squared magnitude of a signal's Fourier transform,
% normalized by the number of frequency samples.
ynoise = fft(xnoise);
ynoiseshift = fftshift(ynoise);    
power = abs(ynoiseshift).^2/n; 

end