function [pacmat, freqvec_ph, freqvec_amp, pac_angles, comodulograms, modulation_indices] = find_pac_nofilt (sig_pac, Fs,...
    measure, sig_mod, ph_freq_vec, amp_freq_vec, plt, waitbar, calc_comodulograms, width, nfft,...
        dataname, sig_pac_name, sig_mod_name)
%function [pacmat, freqvec_ph, freqvec_amp] = find_pac_nofilt (sig_pac, Fs,...
%   measure, sig_mod, ph_freq_vec, amp_freq_vec, plt, waitbar, width, nfft,...
%       dataname, sig_pac_name, sig_mod_name)
%
% This function calculates a matrix of PAC values using either the ESC, MI 
% or CFC measure. 
% It assumes that the input is a prefiltered signal since this function 
% does not include any filtering. The output is a matrix of PAC values and
% a plot of this matrix (depending on the value of the 'plt' argument).
% This signal can only take single vector signals or matrices (for multiple
% trials)
% 
%
% Basic function call:
% function pacmat = find_pac_nofilt(sig_pac, Fs, measure)
%
% REQUIRED INPUTS:
%  sig_pac - signal suspected of containing PAC
%  Fs - sampling frequency
%  measure - measure to be used - it should be: 'esc', 'mi' or 'cfc'
%
% The function can be executed with many optional inputs:
%  function pacmat = find_pac_nofilt(sig_pac, Fs, measure, ...
%    sig_mod, ph_freq_vec, amp_freq_vec, plt, width, nfft, dataname,...
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
% Edited by David Stanley, March 2015

% Checks of input variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 4
    sig_mod = sig_pac;
end
if nargin < 5
    ph_freq_vec = 1:5:101;
end
if nargin < 6
    amp_freq_vec = 1:5:101;
end
if nargin < 7
    plt = 'y';
end
if nargin < 8 
    waitbar = 0;
end
if nargin < 9
    calc_comodulograms = 0;
end
if nargin < 10
    width = 7;
end
if nargin < 11
    nfft = ceil(Fs/(diff(ph_freq_vec(1:2))));
end
if nargin < 12
    dataname = '';
end
if nargin < 13
    sig_pac_name = '';
end
if nargin < 14
    sig_mod_name = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up some parameters for clarity/storage
xbins = ceil((max(ph_freq_vec) - min(ph_freq_vec))/(diff(ph_freq_vec(1:2))));
ybins = ceil((max(amp_freq_vec) - min(amp_freq_vec))/(diff(amp_freq_vec(1:2))));
cent_freq_vec = zeros(xbins,1);
cent_freq_vec2 = zeros(ybins, 1);


if isa(sig_pac, 'cell')
    % Data represents filtered signals
    num_trials = size(sig_pac{1,1},2); 
else
    % Data represents signals 
    num_trials = size(sig_pac,2);
    
end

if waitbar == 1
    counter = 0;
    countermax = xbins*ybins;
    fprintf('\nCalculating PAC values\n');
end

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

% Calculate PAC measures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmp(measure, 'esc') || strcmp(measure, 'mi'))
pacmat = zeros(ybins, xbins);
pac_angles = zeros(ybins, xbins);
comodulograms = cell(ybins, xbins);
modulation_indices = zeros(ybins, xbins);
    
    for i = 1:ybins
        for j = 1:xbins
            if i == 16 && j == 10
                % Dummy code to allow one to inspect different values via breakpoints
                %keyboard
            end
            % Calculate matrix of PAC values
            if strcmp(measure, 'esc')
                pacmat(i,j) = esc_measure(sig_mod{1,j}, ...
                    sig_pac{1,i}, 'y');
            end
            
            if strcmp(measure, 'mi')
                % Pacmat full of raw mi values, not yet normalized
                % [pacmat(i,j), pac_angles(i,j), comodulograms{i,j}, modulation_indices(i,j)] = mi_measure(sig_mod{1,j}, ...
                %     sig_pac{1,i}, calc_comodulograms);
                [pacmat(i,j), pac_angles(i,j), comodulograms{i,j}, modulation_indices(i,j)] = mi_measure(sig_mod{1,j}, ...
                    sig_pac{1,i}, calc_comodulograms, Fs);
            end
            
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
    
end



if strcmp(measure, 'cfc')
    
    % Calculate matrix of PAC values (these have been averaged over trials)
    [pacmat, freqvec_ph, freqvec_amp] = cfc_measure(sig_pac, ...
             sig_mod, 'y',ph_freq_vec, freqvec_amp, Fs, nfft, width, waitbar);
              
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(plt, 'y')
    pac_plot_fun(pacmat, freqvec_ph, freqvec_amp, measure, sig_pac_name,...
                 sig_mod_name, dataname, sig_pac, Fs, sig_mod);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     
 

