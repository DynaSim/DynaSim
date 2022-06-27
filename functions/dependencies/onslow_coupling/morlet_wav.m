function y = morlet_wav(f,t,width)
% function y = morlet_wav(f,t,width)
% 
% Create a Morlet wavelet 'y' with frequency resolution 'f' and temporal 
% resolution 't'. The wavelet will be normalized so the total energy is 1.
% The 'width' defines the temporal and frequency resolution for the given
% centre frequency 'f' by determining the number of cycles of the wavelet
% itself (see Tallon-Baudry et al., J. Neurosci. 15, 722-734 (1997) or 
% Event-Related Potentials: A Methods Handbook, Handy (editor), MIT Press, 
% (2005))
%
% INPUTS:
% f - frequency to filter at 
% s - signal to be filtered
% Fs - sampling frequency of s
% width - parameter which defines the mother wavelet (>= 5 is suggested)
%
% Author: Ole Jensen, August 1998 


sf = f/width;
st = 1/(2*pi*sf);
A = 1/sqrt((st*sqrt(pi)));
y = A*exp(-t.^2/(2*st^2)).*exp(1i*2*pi*f.*t);

