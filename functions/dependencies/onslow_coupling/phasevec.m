function y = phasevec(f,s,Fs,width)
% function y = phasevec(f,s,Fs,width)
%
% Returns a vector 'y' containing the instantaneous phase values of signal 
% 's' filtered for frequency 'f' (via convolution with a complex Morlet 
% wavelet described by 'width')
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

t=-3.5*st:dt:3.5*st; 
m = morlet_wav(f,t,width);
y = conv(s,m');
y = angle(y);
y = y(ceil(length(m)/2):length(y)-floor(length(m)/2));