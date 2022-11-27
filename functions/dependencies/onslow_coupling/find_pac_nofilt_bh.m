function [pacmat, pac_angles, comodulograms, modulation_indices] = find_pac_nofilt_bh(...
  filtered_ph, filtered_amp, Fs, calc_comodulograms)
% This function calculates a matrix of PAC values using the Modulation
% Index method, but built using the Hilbert transform, NOT wavelets.
%
% NOTE: This is heavily adapted from Angela Onslow's phase-amplitude
% coupling code of 2010. Original comments follow.
%
%% Calculate PAC measures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize
comodulograms = cell(1,1);

%% Hilbert transforms of filtered signals (as opposed to wavelets)
analytic_ph = hilbert(filtered_ph{1,1});
analytic_amp = hilbert(filtered_amp{1,1});

phase_hilbert = angle(analytic_ph);
amp_hilbert = abs(analytic_amp);

[pacmat, pac_angles, comodulograms{1,1}, modulation_indices] = mi_measure(...
  phase_hilbert, amp_hilbert, calc_comodulograms, Fs);
