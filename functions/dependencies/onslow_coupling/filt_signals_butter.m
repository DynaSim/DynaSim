function [ph_filt_signals, amp_filt_signals] = filt_signals_butter(sig_ph, sig_amp,...
    Fs, freq_vec_ph, freq_vec_amp, filter_order_ph, filter_order_amp)
% Returns cell arrays 'ph_filt_signals' and 'amp_filt_signals'. These
% contain 'sig_amp' and 'sig_ph' filtered according to Butterworth filters.
%
% NOTE: This is heavily adapted from Angela Onslow's phase-amplitude coupling code of 2010.
%
% INPUTS:
% sig_amp - signal which will be analysed as containing the higher frequency,
% modulated PAC signal
% sig_ph - signal which will be analysed as containing the lower frequency,
% modulating signal
% Fs - sampling frequency, in Hz
% freq_vec_ph - vector of frequencies to filter sig_ph for, in Hz
% freq_vec_amp - vector of frequencies to filter sig_amp for, in Hz
% measure - PAC measure which will be calculated using these filtered
% signals
% width - number of cycles which will define the mother Morlet wavelet
%
% OUTPUTS:
% ph_filt_signals - 1 x (number of bins determined by freq_vec_ph) cell 
% array, each element has the same dimensions as the original signals
% amp_filt_signals - 1 x (number of bins determined by freq_vec_amp) cell 
% array, each element has the same dimensions as the original signals
%
% Author: Angela Onslow, May 2010


total_num_dp = size(sig_amp,1);
num_trials = 1:size(sig_amp,2);

ph_freq_low = min(freq_vec_ph);
ph_freq_high = max(freq_vec_ph);
amp_freq_low = min(freq_vec_amp);
amp_freq_high = max(freq_vec_amp);

% Create structures to store filtered signals
ph_filt_signals = cell(1,1);
amp_filt_signals = cell(1,1);

ph_filt_signals{1,1} = zeros(total_num_dp,max(num_trials));
amp_filt_signals{1,1} = zeros(total_num_dp,max(num_trials));

% Filter and store filtered signals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i3=1:max(num_trials)
    ph_filt_signals{1,1}(:,i3) =  apply_filter_bh(ph_freq_low, ph_freq_high, sig_ph(:,i3),...
                                                  Fs, filter_order_ph);
    amp_filt_signals{1,1}(:,i3) = apply_filter_bh(amp_freq_low, amp_freq_high, sig_amp(:,i3),...
                                                  Fs, filter_order_amp);
end
