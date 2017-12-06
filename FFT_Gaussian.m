% Convert a Gaussian pulse from the time domain to the frequency domain.

Fs  = 100;              % Sampling frequency
t   = -0.5:1/Fs:0.5;    % Time vector 
L   = length(t);        % Signal length
X   = 1/(4*sqrt(2*pi*0.01))*(exp(-t.^2/(2*0.01)));

if 1                    % Plot the pulse in the time domain
    plot(t,X)
    title('Gaussian Pulse in Time Domain')
    xlabel('Time (t)')
    ylabel('X(t)')
end

% To use the fft function to convert the signal to the frequency domain, 
% first identify a new input length that is the next power of 2 from the 
% original signal length. This will pad the signal X with trailing zeros 
% in order to improve the performance of fft.
n = 2^nextpow2(L);
Y = fft(X,n);           % Y = fft(X,n) returns the n-point discrete Fourier transform  (DFT)
f = Fs*(0:(n/2))/n;
P = abs(Y/n);

plot(f,P(1:n/2+1)) 
title('Gaussian Pulse in Frequency Domain')
xlabel('Frequency (f)')
ylabel('|P(f)|')