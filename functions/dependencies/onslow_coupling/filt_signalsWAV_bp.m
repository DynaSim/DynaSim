function [bp_sig1, bp_sig2] = filt_signalsWAV_bp(sig1,sig2,...
    Fs, freq_vec, width)
% function [bp_sig1, bp_sig2] = filt_signalsWAV_bp(sig1,sig2,...
%    Fs, freq_vec, width)
%
% This function band-pass filters sig1 and sig2 using a wavelet filter.
% Bandwidth of filter assumed equal to the difference between max(freq_vec)
% and min(freq_vec)
%
% Inputs: 
%
% sig1 - first signal to be band-pass filtered
%
% sig2 - second signal to be band-pass filtered
%
% Fs - sampling frequency, in Hz
%
% freq_vec - vector describing the band-pass region, e.g to filter between
% 5 and 10 Hz freq_vec should be equal to [5 10]
%
% width - parameter defining the mother wavelet used for filtering, 7 is a
% good default value
%
% Author: Angela Onslow, 2010


num_trials = size(sig1,2);

freq_low = min(freq_vec);
freq_high = max(freq_vec);

% Create structures to store filtered signals
% bp_sig1 = zeros(1,num_trials);
% bp_sig2 = zeros(1,num_trials);
bp_sig1 = [];
bp_sig2 = [];


% Filter and store filtered signals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bandpass filter both signals

for i=1:num_trials
    bp_sig1(:,i) = bpvec((freq_low + floor((freq_high - freq_low)/2)), sig1(:,i), Fs, width);
    bp_sig2(:,i) = bpvec((freq_low + floor((freq_high - freq_low)/2)), sig2(:,i), Fs, width);
end


