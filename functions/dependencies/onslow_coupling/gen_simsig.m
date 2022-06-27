function[simsig simsigmod snr] = gen_simsig(sigy1_freq, sigy2_freq, Fs, ...
    numsec, noise_ratio, amp_ratio)
% function[simsig simsigmod snr] = gen_simsig(sigy1_freq, sigy2_freq, Fs,...
%   numsec, noise_ratio, amp_ratio)
%
% This function generates a PAC containing signal 'simsig' and a
% corresponding modulating signal 'simsigmod'.
%
% INPUTS:
% sigy1_freq - lower, modulating signal frequency, in Hz
% sigy2_freq - higher, modulated signal frequency, in Hz
% Fs - sampling frequency, in Hz
% numsec - length of generated signals, in seconds
% noise_ratio - standard deviation of noise to be added to generated
% signals
% amp_ratio - scaling factor with which to vary the amplitude of the higher
% frequency signal
%
% OUTPUTS:
% simsig - PAC containing signal, plus noise
% simsigmod - modulating frequency signal, plus noise
% snr - value for the signal-to-noise ratio
%
% Author: Angela Onslow, May 2010

if sigy2_freq > Fs/2
    printf('Warning, your signal does not meet the Nyquist criterion: the highest frequency present should be < Fs/2'); 
end

sigx = 1:(Fs*numsec); 

% Lower frequency modulator signal
sigy1 = sin((sigy1_freq/Fs)*2*pi*sigx);

% Higher frequency modulated signal
sigy2 = (sin((sigy2_freq/Fs)*2*pi*sigx)).*amp_ratio;

% Final PAC signal
sigy3 = (sigy2.*(sigy1+1));

% Create noise signals
noise_sig1 = noise_ratio.*randn(max(sigx),1);
noise_sig2 = noise_ratio.*randn(max(sigx),1);

% Signal containing PAC + noise
simsig = sigy3' + noise_sig2;

% Signal containing low frequency + noise
simsigmod = sigy1' + noise_sig1;

% Calculate SNR
simsigpow = rms_val(sigy3);
noisepow = rms_val(noise_sig2);
snr = 20*log10((simsigpow/noisepow));