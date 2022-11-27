function filtered_signal = apply_filter_bh(bottom_freq, top_freq, signal, Fs, order)
% function filtered_signal = apply_filter_bh(bottom_freq, top_freq, signal, Fs, order)
%
% Returns a vector 'filtered_signal', taken from having applied a
% bandpass Butterworth filter to 'signal'. The filter is built from
% 'bottom_freq'/Nyquist, 'top_freq'/Nyquist, and of order 'order'.
%
% NOTE: This is heavily adapted from Angela Onslow's phase-amplitude
% coupling code of 2010. Original comments follow.
%
% NOTE: This function is a modification of code written Ole Jensen, August
% 1998, see ENERGYVEC
%
% Author: Angela Onslow, May 2010

nyq = Fs/2;

% Signal has already been detrended before being sent to the "find_pac..." function
% detrended_signal = detrend(signal)

[coeffs_b, coeffs_a] = butter(order, [bottom_freq/nyq, top_freq/nyq], 'bandpass');

filtered_signal = filtfilt(coeffs_b, coeffs_a, double(signal));
