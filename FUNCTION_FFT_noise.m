function [fshift,power] = FUNCTION_FFT_noise(xnoise, delta_t)

Fs      = 1/delta_t;
n       = 2^nextpow2(length(xnoise));
nn      = n*10;
fshift  = Fs*(0:(n/2))/nn;
% fshift  = transpose(fshift);

% Power is the squared magnitude of a signal's Fourier transform,
% normalized by the number of frequency samples.
ynoise = fft(xnoise, nn);
ynoiseshift = fftshift(ynoise);    
power = abs(ynoiseshift((nn/2+1):(nn/2+n/2+1))).^2/nn; 

end

% synchronize the file to the data folder 
% rsync /Volumes/DataSSD/MATLAB_codes/Project180129-FT/FUNCTION_FFT_noise.m /Volumes/DataSSD/OneDrive\ -\ UNSW/Hermite_Decomposition/ESO_HARPS/code