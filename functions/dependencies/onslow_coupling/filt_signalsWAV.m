function [ph_filt_signals, amp_filt_signals] = filt_signalsWAV(sig1,sig2,...
    Fs, ph_freq_vec, amp_freq_vec, measure, width)
% function [ph_filt_signals, amp_filt_signals] = filt_signalsWAV(sig1,sig2,...
% Fs, ph_freq_vec, amp_freq_vec, measure, width)
%
% Returns cell arrays 'ph_filt_signals' and 'amp_filt_signals'. These 
% contain 'sig1' and 'sig2' correctly filtered in order to calculate the
% PAC measure, either ESC or MI given by 'measure'. The frequencies
% to filter at are determined by the vectors 'ph_freq_vec' and
% 'amp_freq_vec'. The frequency range is defined by the minimum and the
% maximum value and the interval, or bandwith, between values. For example
% if 'ph_freq_vec' is defined as 1:5:101 ([1, 6, 11, 16, 21...101]) then
% the centre frequencies that will be filtered for are [3, 8, 13, 18...98]
% i.e. the centre of the 5 Hz bins from 1 Hz to 101 Hz. Filtering is
% achieved via convolution with a complex Morlet wavelet (number of cycles 
% defined by 'width' parameter). 
%
% INPUTS:
% sig1 - signal which will be analysed as containing the higher frequency,
% modulated PAC signal
% sig2 - signal which will be analysed as containing the lower frequency,
% modulating signal
% Fs - sampling frequency, in Hz
% ph_freq_vec - vector of frequencies to filter sig2 for, in Hz
% amp_freq_vec - vector of frequencies to filter sig1 for, in Hz
% measure - PAC measure which will be calculated using these filtered
% signals
% width - number of cycles which will define the mother Morlet wavelet
%
% OUTPUTS:
% ph_filt_signals - 1 x (number of bins determined by ph_freq_vec) cell 
% array, each element has the same dimensions as the original signals
% amp_filt_signals - 1 x (number of bins determined by amp_freq_vec) cell 
% array, each element has the same dimensions as the original signals
%
% Author: Angela Onslow, May 2010

total_num_dp = size(sig1,1);
num_trials = 1:size(sig1,2);

phfreq_low = min(ph_freq_vec);
phfreq_high = max(ph_freq_vec);
phfreq_bw = diff(ph_freq_vec(1:2));
ampfreq_low = min(amp_freq_vec);
ampfreq_high = max(amp_freq_vec);
ampfreq_bw = diff(amp_freq_vec(1:2));

xbins = ceil((phfreq_high - phfreq_low)/phfreq_bw);
ybins = ceil((ampfreq_high - ampfreq_low)/ampfreq_bw);


% Create structures to store filtered signals
ph_filt_signals = cell(1,xbins);
amp_filt_signals = cell(1,ybins);

for i2=1:xbins
    ph_filt_signals{1,i2} = zeros(total_num_dp,max(num_trials));
end


for i2=1:ybins
    amp_filt_signals{1,i2} = zeros(total_num_dp,max(num_trials));
end

% Filter and store filtered signals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filter sig2 for phase freq
for i2=1:xbins
    
    upper_bin = phfreq_low+(i2*phfreq_bw);
    lower_bin = upper_bin-phfreq_bw;
     
    if strcmp(measure, 'esc')
        for i3=1:max(num_trials)
            ph_filt_signals{1,i2}(:,i3) = bpvec((lower_bin + floor((upper_bin- lower_bin)/2)), sig2(:,i3), Fs, width);
        end
    else
        for i3=1:max(num_trials)
            ph_filt_signals{1,i2}(:,i3) = phasevec((lower_bin + floor((upper_bin- lower_bin)/2)),sig2(:,i3), Fs,width);
        end
    end
    
       
end

% Filter sig1 for amplitude freq
for i2=1:ybins
    
    upper_bin = ampfreq_low+(i2*ampfreq_bw);
    lower_bin = upper_bin-ampfreq_bw;
            
    for i3=1:max(num_trials)
        amp_filt_signals{1,i2}(:,i3) = ampvec((lower_bin + floor((upper_bin- lower_bin)/2)),sig1(:,i3), Fs, width);

    end
    
  
end
