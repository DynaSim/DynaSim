function [pacmat, freqvec_ph, freqvec_amp] = find_pac_nofilt_win(sig_pac, Fs,...
    measure, win_length, overlap, sig_mod, ph_freq_vec, amp_freq_vec, plt, waitbar, width, nfft,...
        dataname, sig_pac_name, sig_mod_name)
%function [pacmat, freqvec_ph, freqvec_amp] = find_pac_nofilt_win(sig_pac, Fs,...
%   measure, win_length, overlap, sig_mod, ph_freq_vec, amp_freq_vec, plt, waitbar, width, nfft,...
%       dataname, sig_pac_name, sig_mod_name)
%
% This function calculates a matrix of PAC values using either the ESC, MI 
% or CFC measure: the x-axis shows the time evolution of PAC, the y-axis
% shows the higher, modulated frequency (the lower, modulating frequency is
% filtered in one particular band only).
% It assumes that the input is a prefiltered signal since this function 
% does not include any filtering. The output is a matrix of PAC values and
% a plot of this matrix (depending on the value of the 'plt' argument).
%
% Basic function call:
% function pacmat = find_pac_nofilt(sig_pac, Fs, measure, win_length, overlap)
%
% REQUIRED INPUTS:
%  sig_pac - signal suspected of containing PAC
%  Fs - sampling frequency
%  measure - measure to be used - it should be: 'esc', 'mi' or 'cfc'
%  win_length - length of time windows, in seconds
%  overlap - amount of overlap of time windows, in seconds
%
% The function can be executed with many optional inputs:
%  function pacmat = find_pac_nofilt(sig_pac, Fs, measure, win_length,
%  overlap, sig_mod, ph_freq_vec, amp_freq_vec, plt, width, nfft, dataname,...
%       sig_pac_name, sig_mod_name)
%
% OPTIONAL INPUTS:
%  sig_mod - signal containing modulating frequency; default = sig_pac
%  ph_freq_vec - range of frequencies for modulating signal; default = 1:5:101
%  amp_freq_vec - range of frequencies for PAC signal; default = 1:5:101
%  plt - flag indicating whether the output should be plotted - it should
%        be 'y' or 'n'; default = 'y'
%  width - width of the wavelet filter; default = 7
%  nfft - the number of points in fft; default = 200
%  dataname - the name of the dataset to be included as a graph title;
%             default = ''
%  sig_pac_name - the name of sig_pac to be printed on the y-axis; default = ''
%  sig_mod_name - the name of sig_mod to be printed on the x-axis; default = ''
%
% Author: Angela Onslow, May 2010 

% Checks of input variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 6
    sig_mod = sig_pac;
end
if nargin < 7
    ph_freq_vec = 1:5:101;
end
if nargin < 8
    amp_freq_vec = 1:5:101;
end
if nargin < 9
    plt = 'y';
end
if nargin < 10 
    waitbar = 0;
end
if nargin < 11
    width = 7;
end
if nargin < 12
    nfft = Fs;
end
if nargin < 13
    dataname = '';
end
if nargin < 14
    sig_pac_name = '';
end
if nargin < 15
    sig_mod_name = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up some parameters for clarity/storage
ybins = ceil((max(amp_freq_vec) - min(amp_freq_vec))/(diff(amp_freq_vec(1:2))));


% % lf_bins = ceil((max(ph_freq_vec) - min(ph_freq_vec))/(diff(ph_freq_vec(1:2))));
% % cent_freq_vec = zeros(xbins,1);
% % cent_freq_vec2 = zeros(ybins, 1);


if isa(sig_pac, 'cell')
    % Data represents filtered signals
    num_trials = size(sig_pac{1,1},2);
    %xbins = (size(sig_pac{1,1},1))/ceil(win_length*Fs);
else
    % Data represents signals - CRC measure
    num_trials = size(sig_pac,2);
    %xbins = (size(sig_pac,1))/ceil(win_length*Fs);
    
end

if waitbar == 1
    counter = 0;
    countermax = ybins;
    fprintf('\nCalculating PAC values\n');
end

%if strcmp(measure,'esc')||strcmp(measure, 'mi')
    xbins = ceil((max(ph_freq_vec) - min(ph_freq_vec))/(diff(ph_freq_vec(1:2))));
    
for i =1:xbins
    upper_bin = min(ph_freq_vec)+i*(diff(ph_freq_vec(1:2)));
    lower_bin = upper_bin-(diff(ph_freq_vec(1:2)));
    cent_freq_vec(i) =  lower_bin + floor((upper_bin- lower_bin)/2);
end


for i =1:ybins
    upper_bin = min(amp_freq_vec)+i*(diff(amp_freq_vec(1:2)));
    lower_bin = upper_bin-(diff(amp_freq_vec(1:2)));
    cent_freq_vec2(i) =  lower_bin + floor((upper_bin- lower_bin)/2);
end


freqvec_amp = cent_freq_vec2;
freqvec_ph = cent_freq_vec;

%end


% Calculate PAC measures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmp(measure, 'esc') || strcmp(measure, 'mi'))
    %pacmat = zeros(ybins, xbins);
    pacmat = [];
    
    for i = 1:ybins
        
        pacmat(i,:) = window_data(sig_pac{1,i}, Fs, win_length, overlap,...
            measure, sig_mod{1,1});
        
        
        % Display current computational step to user
        if waitbar ==1
            counter = counter+1;
            if counter == 1
                fprintf('%03i%% ', floor((counter/countermax)*100));
            else
                fprintf('\b\b\b\b\b%03i%% ', floor((counter/countermax)*100));
            end
            if counter == countermax
                fprintf('\n');
            end
        end
    end
end



if strcmp(measure, 'cfc')
    
    [pacmat, freqvec_ph, freqvec_amp] = window_data(sig_pac,Fs, ...
    win_length,overlap, measure, sig_mod,ph_freq_vec, freqvec_amp, nfft, width, waitbar);
              
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(plt, 'y')
    pac_plot_fun(pacmat, freqvec_ph, freqvec_amp, measure, sig_pac_name,...
                 sig_mod_name, dataname, sig_pac, Fs, sig_mod);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     
 

