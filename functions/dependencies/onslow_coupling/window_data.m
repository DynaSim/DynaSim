function [out_mat, freqvec_ph, freqvec_amp] = window_data(sig,Fs, ...
    win_length,overlap, measure, sig2,ph_freq_vec, freqvec_amp, nfft, width, waitbar)
% function [out_mat, freqvec_ph, freqvec_amp] = window_data(sig,Fs, ...
%    win_length,overlap, measure, sig2,ph_freq_vec, freqvec_amp, nfft, width, waitbar)
%
% This function windows the signals passed as 'sig' and (optionally) 'sig2' and then calculates 
% the PAC measure of choice on each windowed segment. 'sig' is assumed to contained PAC whilst 
% 'sig2' is assumed to contain the lower, modulating frequency responsible. Trials must be 
% passed as columns of a matrix. In order to plot the time-varying evolution of PAC on the
% x-axis the input argument 'ph_freq_vec' should be 1 number = to the lower, modulating
% frequency of interest.
%
% Inputs:
%
% sig - signal believed to contain PAC
%
% Fs - sampling frequency
%
% win_length - length of desired windows in seconds (N.B. if your win_length = 0.5s then the 
% lowest frequency you can accurate measure in your data is 2 Hz i.e. 1/win_length)
%
% overlap - length of window overlap in seconds
%
% measure - PAC measure of choice: 'ESC', 'MI' or 'CFC
%
% sig2 - signal believed to contained the lower, modulating frequency 
%
% ph_freq_amp - lower, modulating frequency of interest
%
% freqvec_amp - range of frequencies to look at in PAC signal
%  
% nfft - the number of points in fft
%
% width - width of the wavelet filter
%
% waitbar - display progress in the command window
%
% Outputs:
%
% out_mat - vector of PAC values, each columns corresponds to the PAC value for one time window
%
% freqvec_ph - vector of values at which the lower, modulating frequency has been examined 
% (mainly as a check)
%
% freqvec_amp - vector of values at which the higher, modulated frequency has been examined
%
% Author: Angela Onslow, April 2011



% sig assumed sig_pac
% sig2 assumed sig_mod

if nargin < 5
    sig2 = sig;
end

if (strcmp(measure, 'cfc')) && (nargin < 11)
    sprintf('function window_data requires additional arguments when measure = cfc');
    return
end

% This function assumes trials are passed as columns, i.e sig contains more
% rows than columns, if this is not the case the signal will be
% automatically transposed
if size(sig,1)<size(sig,2)
    sig = sig';
end

if size(sig2,1)<size(sig2,2)
    sig2 = sig2';
end

if (size(sig,2) ~= size(sig2,2))
    sprintf('Error - Signals must have the same number of trials')
    return
end

%# compute indices
nx = size(sig,1);
ntrials = size(sig,2);
nwind = ceil(win_length*Fs);
noverlap = ceil(overlap*Fs);
idx = bsxfun(@plus, (1:nwind)', 1+(0:(fix((nx-noverlap)/(nwind-noverlap))-1))*(nwind-noverlap))-1;


if strcmp(measure,'cfc')
    out_mat = [];
else
    out_mat = zeros(1,size(idx,2));
end

%# loop over sliding windows
for k=1:size(idx,2)
    
    win_mat = [];
    win_mat2 = [];
    
    %# loop over trials
    for j = 1:ntrials
        
        slidingWindow = sig(idx(:,k),j);
        slidingWindow2 = sig2(idx(:,k),j);
        
        win_mat = [win_mat, slidingWindow];
        win_mat2 = [win_mat2, slidingWindow2];
        
    end
    
    if strcmp(measure, 'esc')
        
        out_mat(1,k) = esc_measure(win_mat2,win_mat, 'y');
    
    end
    
    if strcmp(measure, 'mi')
        
        out_mat(1,k) = mi_measure(win_mat2,win_mat);
  
    end
    
    if strcmp(measure, 'cfc')
        
        [out_mat(:,k), freqvec_ph, freqvec_amp] = cfc_measure(win_mat,win_mat2,...
            'y',ph_freq_vec, freqvec_amp,Fs, nfft, width, 0);
        
        
         if waitbar ==1
            if k ==1
                fprintf('%03i%% ', floor(k/size(idx,2)*100));
            else
                fprintf('\b\b\b\b\b%03i%% ', floor(k/size(idx,2)*100));
            end
            if k == size(idx,2)
                fprintf('\n');
            end
        end
   
    end
end