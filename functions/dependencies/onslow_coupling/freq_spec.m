function [spec maxfreq specpeak f] =  freq_spec(sig, Fs, plt, x_max)
% function [spec maxfreq specpeak] =  freq_spec(sig, Fs, plt, x_max)
%
% This function plots the frequency spectrum (aka power spectrum)of a 
% signal and returns the frequency with the largest power as 'maxfreq'
%
% INPUTS:
% sig - signal to be investigated
%
% Fs - sampling frequency, in Hz
%
% plt - passed as a string, if 'y' then the frequency spectrum will be 
% plotted
%
% x_max - maximum frequency to plot on x-axis of spectrum
%
% OUTPUTS:
% spec - the frequency spectrum values passed as a row vector, 1st element
% is the power of the signal components of 1 Hz frequency and so on
%
% maxfreq - the frequency value with the largest power contribution
%
% specpeak - the amplitude of the predominant oscillation
%
% Author: Angela Onslow, December 2009

if nargin < 4
    x_max = floor(Fs/2);
end

L = length(sig);

% NFFT - number of points describing the fft calculation which will
% produce a desired frequency resolution, the fft algorithm can run much
% faster if this number is a power of 2.
% Frequency Resolution = Sample Rate / FFT Points (NFFT) = 
% Frequency range / Number of bins
NFFT = 2^nextpow2(L);

% Vector of frequencies at which the FFT will estimate power.
% The maximum frequency which can be established from a sampled signal is 
% Fs/2 Hz due to the Nyquist theorem
f = Fs/2*linspace(0,1,NFFT/2+1);
if max(f) > x_max
    [C,ia] = intersect(round(f),x_max);
    f = f(1:ia);
end
    

% Calculate FFT, normalize such that area under the curve = 1
signalFFT = fft(sig, NFFT)/L;

% Signal spectrum is equal to the first half of signalFFT, the second half
% contains information at negative frequencies (a mathematical quirk, this
% information is redundant), then take the absolute value (real part) of
% signalFFT since these values are complex numbers and square to get power
% (Power = Energy/Time, signal energy equals signal amplitude squared)
spec = 2*abs(signalFFT(1:NFFT/2+1));
if length(spec) > length(f)
    spec = spec(1:ia);
end

if nargin >2
    
if strcmp(plt, 'y')
    
% Plot spectrum    
figure
plot(f, spec);
title('Frequency Spectrum')
xlabel('Frequency Hz')
ylabel('Amplitude')
%semilogy(f, spec);
%loglog(f,spec);
xlim([0 x_max])

end
end

% Find maximum frequency component and return as 'maxfreq'
[specpeak ind] = max(spec);
maxfreq = f(1, ind);
