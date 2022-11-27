function y = ampvec_bh(bottom_freq, top_freq, s, Fs, order)
% function y = ampvec(f,s,Fs,width)
%
% Returns a vector 'y' containing the instantaneous amplitude values of 
% signal 's' filtered for frequency 'f' (via convolution with a complex 
% Morlet wavelet described by 'width')
%
% INPUTS:
% f - frequency to filter at 
% s - signal to be filtered
% Fs - sampling frequency of s
% width - parameter which defines the mother wavelet (this is then 
% scaled and translated in order to filter for different frequencies, 
% >= 5 is suggested, see Tallon-Baudry et al., J. Neurosci. 15, 722-734 
% (1997))
%
% NOTE: This function is a modification of code written Ole Jensen, August
% 1998, see ENERGYVEC
%
% Author: Angela Onslow, May 2010

dt = 1/Fs;
sf = f/width;
st = 1/(2*pi*sf); 

t=-(width/2)*st:dt:(width/2)*st; 
m = morlet_wav(f,t,width);
y = conv(s,m');
y = abs(y);
y = y(ceil(length(m)/2):length(y)-floor(length(m)/2));



[slow_b, slow_a] = butter(2, [0.5/nyq, 2/nyq], 'bandpass');
% [slow_b, slow_a] = butter(2, [0.5/nyq, 2/nyq], 'bandpass');
% [slow_b, slow_a] = butter(2, [0.5/nyq, 1.5/nyq], 'bandpass');
% [slow_b, slow_a] = butter(2, [1.5/nyq, 2/nyq], 'bandpass');

[alpha_b, alpha_a] = butter(2, [8/nyq, 12/nyq], 'bandpass');
% [alpha_b, alpha_a] = butter(2, [8/nyq, 14/nyq], 'bandpass');

% figure(nn)
% freqz(slow_b, slow_a, [], fs)
% % xlim([0 50])

filtered_slow = filtfilt(slow_b, slow_a, double(detrended_lfp));
filtered_alpha = filtfilt(alpha_b, alpha_a, double(detrended_lfp));

