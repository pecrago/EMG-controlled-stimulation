function [noiseOut] = signalnoise(length, f_sample, low, high)
% noiseOut is a column vector of (length,1) of random noise filtered by a 2nd order
% butterworth bandpass between the low and high cutoffs
% f_sample is the sample rate
% low, high, and f_sample are all in Hz, or all in rad/s
n = 2; % filter order
Wn = [low high]/f_sample/2; % 1 is the Nyquist rate = half the sample rate
ftype = 'bandpass';
[b,a] = butter(n, Wn, ftype);
noiseIn = randn(length, 1); % generate a random signal
noiseOut = filter(b, a, noiseIn); % filter the signal
