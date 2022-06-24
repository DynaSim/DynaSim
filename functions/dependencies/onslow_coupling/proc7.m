function[simsig simsigmod snr] = proc7(sigy1_freq, sigy2_freq, Fs, numsec, noise_ratio, amp_ratio, phshift)
%function[simsig simsigmod snr] = proc7(sigy1_freq, sigy2_freq, Fs, numsec, noise_ratio, amp_ratio, phshift)
%
%This function creates 2 signals, one (simsig) which contains a PAC wave plus the
%lower modulating frequency (this second component can be taken out easily within
%the code but is currently not available as an argument) and one
%(simsigmod) which contains only the lower modulating frequency. Both
%signals have some level of noise.
%
%This code differs from proc5 in allowing the user to specify a phase shift
%which is added to both signals.
%
%Inputs:
%
%sigy1_freq: frequency in Hz of the lower, modulating frequency, e.g. 4
%
%sigy2_freq: frequency in Hz of the higher, modulated frequency, e.g 60
%
%Fs: assumed sampling frequency of the created signals
%
%numsec: number of seconds worth of data you wish to create
%
%noise_ratio: the level of noise which will be added to each signal, e.g.
%0.5
%
%amp_ratio: the ratio of amplitudes between the high and low frequency
%components i.e. between the PAC signal and the modulating signal
%
%phshift: the phase shift to introduce into the signals
%
%Outputs:
%
%simsig: created signal, containing a PAC component, a lower, modulating 
%frequency component + noise
%
%simsigmod: created signal, containing the lower, modulating frequency 
%component + noise
%
%Example usage:
%
%[simsig simsigmod] = proc7(4,60,1017,10,0.5,1,(pi/2))
%
%Author: Angela Onslow, December 2009


sigx = 1:(Fs*numsec); 

%carrier (modulator)
sigy1 = sin((sigy1_freq/Fs)*2*pi*sigx + phshift);

%message (to be modulated)
sigy2 = (sin((sigy2_freq/Fs)*2*pi*sigx)).*amp_ratio;

%final modulated signal
sigy3 = (sigy2.*(sigy1+1));

noise_sig1 = noise_ratio.*randn(max(sigx),1);
noise_sig2 = noise_ratio.*randn(max(sigx),1);

for i=1:1
    
    sigHPCChoice(:,i) = sigy1' + noise_sig1;
    sigHPCForced(:,i) = randn(max(sigx),1);
    sigPFCChoice(:,i) = sigy3' + noise_sig2;
    sigPFCForced(:,i) = randn(max(sigx),1);
    
end

%simsig = (sigHPCChoice(:,1) + amp_ratio*sigPFCChoice(:,1));
%Uncomment the code below in order to remove the lower frequency component 
simsig = sigPFCChoice(:,1);
simsigmod = (sigHPCChoice(:,1));

simsigpow = rms_val(sigy3);
noisepow = rms_val(noise_sig2);
snr = 20*log10((simsigpow/noisepow));

% simsigpow = rms_sq_val(sigy3);
% snr = simsigpow/((noise_ratio)^2);


