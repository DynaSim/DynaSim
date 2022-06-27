function [freq_mat1, freq_mat2, maxfreq_store1, maxfreq_store2] = find_peak_freq_tf (sig_struct1, ...
    sig_struct2, Fs, win_length, overlap, ylim, filt_freq, plt, waitbar, width,...
    dataname, struct1_name, struct2_name)
% This function calculates the peak frequency of two signals over the course
% of several time windows (which can be made to overlap). It plots this
% time course of changes in peak frequency in a variety of different ways in
% order to aid comparison between the signals:
% 
% 1) Heat plot spectrum, two separate plots
% 2) Heat plot difference in two spectrums, single plot
% 3) Line plot maxfreq variation, two separate plots
% 4) Line plot difference in maxfreq, single plot
%
% The two signals may be filtered for a particular frequency band first in
% order to look at changes in peak frequency within a certain band.
%
% Basic function call:
% find_peak_freq_tf (sig_struct1,sig_struct2, Fs, win_length, overlap, ylim)
%
% REQUIRED INPUTS:
%  sig_struct1 - 1st signal under analysis
%  sig_struct2 - 2nd signal under analysis
%  Fs - sampling frequency
%  win_length - length of time windows, in seconds
%  overlap - amount of overlap of time windows, in seconds
%  y_lim - max frequency component to analyse (= limit of y-axis on plots)
%
% The function can be executed with many optional inputs
%  [freq_mat1, freq_mat2, maxfreq_store1, maxfreq_store2] = find_peak_freq_tf (sig_struct1, ...
%    sig_struct2, Fs, win_length, overlap, ylim, filt_freq, plt, waitbar, width,...
%    dataname, struct1_name, struct2_name)
%
% OPTIONAL INPUTS:
%  filt_freq - vector describing a band-pass filter which will be applied
%  to both signals before looking at changes in peak frequency, e.g [5 10]
%  describes a band-pass filtering centred on 7 Hz and encompassing
%  frequencies between 5 and 10 Hz
%  plt - flag indicating whether the output should be plotted - 0 for not plot, 
%  1 for plot option 1) (above) (this is the default setting), 2 for plot 
%  option 2) and so on
%  waitbar - display progress in the command window; 1 = yes, default = 0 
%  (no)
%  width - width of the wavelet filter; default = 7
%  dataname - the name of the dataset to be included as a graph title;
%             default = ''
%  struct1_name - the name of sig_struct1 to be printed on the y-axis; default = ''
%  struct2_name - the name of sig_struct2 to be printed on the x-axis; default = ''
%
% OUTPUT:
%  freq_mat1 - matrix of values of power at each frequency (rows) for each
%  time window (columns) for sig_struct1, therefore one column = power spectrum
%  for that time window
%  freq_mat2 - matrix of values of power at each frequency (rows) for each
%  time window (columns) for sig_struct2
%  maxfreq_store1 - vector of values of peak frequency for each time window
%  for sig_struct1
%  maxfreq_store2 - vector of values of peak frequency for each time window
%  for sig_struct2
%
% Author: Angela Onslow, May 2011

% Checks of input variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 7
    filt_freq = 0;
end
if nargin < 8
    plt = 1;
end
if nargin < 9
    waitbar = 0;
end
if nargin < 10
    width = 7;
end
if nargin < 11
    dataname = '';
end
if nargin < 12
    struct1_name = '';
end
if nargin < 13
    struct2_name = '';
end

% Check data is columnwise
if size(sig_struct1,1)<size(sig_struct1,2)
    sig_struct1 = sig_struct1';
end

if size(sig_struct2,1)<size(sig_struct2,2)
    sig_struct2 = sig_struct2';
end

if (size(sig_struct1,2) ~= size(sig_struct2,2))
    sprintf('Error - Signals must have the same number of trials')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up global params here
num_trials = size(sig_struct1,2);
nx = size(sig_struct1,1);
nwind = ceil(win_length*Fs);
noverlap = ceil(overlap*Fs);
cut_sec = mod(nx,nwind);
sig_struct1 = sig_struct1(1:(nx-cut_sec),:);
sig_struct2 = sig_struct2(1:(nx-cut_sec),:);
% New nx
nx = size(sig_struct1,1);
idx = bsxfun(@plus, (1:nwind)', 1+(0:(fix((nx-noverlap)/(nwind-noverlap))-1))*(nwind-noverlap))-1;

if num_trials == 1
    freq_mat1=[];
    freq_mat2=[];
    maxfreq_store1 = [];
    maxfreq_store2 = [];
else
    freq_mat1cell = cell(((2^nextpow2(nwind))/2)+1,size(idx,2));
    freq_mat2cell = cell(((2^nextpow2(nwind))/2)+1,size(idx,2));
    
end

% Filter signals if necessary
if filt_freq ~= 0
    
    [sig_struct1, sig_struct2] = filt_signalsWAV_bp(sig_struct1,sig_struct2,...
        Fs, filt_freq, width);
end

% Window data

if num_trials == 1
    if waitbar == 1
        fprintf('\nProgress:\n');
    end
        
    
    
    for k=1:size(idx,2)
        
        
        
        
        slidingWindow1 = sig_struct1(idx(:,k),1).*hamming(nwind,'periodic');
        slidingWindow2 = sig_struct2(idx(:,k),1).*hamming(nwind,'periodic');
        
        % Calculate spectrum and store
        
        [spec1 maxfreq1 specpeak1 f1] =  freq_spec(slidingWindow1, Fs,'n', ylim);
        [spec2 maxfreq2 specpeak2 f2] =  freq_spec(slidingWindow2, Fs,'n', ylim);
        
        
        freq_mat1(:,k) = spec1;
        freq_mat2(:,k) = spec2;
        
        maxfreq_store1(1,k) = maxfreq1;
        maxfreq_store2(1,k) = maxfreq2;
        
        % Display current computational step to user
        if waitbar == 1
            if k == 1
                fprintf('%03i%% ', floor((k/(size(idx,2)))*100));
            else
                fprintf('\b\b\b\b\b%03i%% ', floor((k/(size(idx,2)))*100));
            end
            if k == size(idx,2)
                fprintf('\n');
            end
        end
    end
    
else
    for c = 1:num_trials
        if waitbar == 1
            fprintf('\nTrial_%i:\n', c);
        end
        for k=1:size(idx,2)
            slidingWindow1 = sig_struct1(idx(:,k),c).*hamming(nwind,'periodic');
            slidingWindow2 = sig_struct2(idx(:,k),c).*hamming(nwind,'periodic');
            
            [spec1 maxfreq1 specpeak1 f1] =  freq_spec(slidingWindow1, Fs,'n', ylim);
            [spec2 maxfreq2 specpeak2 f2] =  freq_spec(slidingWindow2, Fs,'n', ylim);
            
            freq_mat1(:,k) = spec1;
            freq_mat2(:,k) = spec2;
            
            maxfreq_store1(c,k) = maxfreq1;
            maxfreq_store2(c,k) = maxfreq2;
            %maxfreq_store3(c,k) = abs(maxfreq1-maxfreq2);
            
            % Display current computational step to user
        if waitbar == 1
            if k == 1
                fprintf('%03i%% ', floor((k/(size(idx,2)))*100));
            else
                fprintf('\b\b\b\b\b%03i%% ', floor((k/(size(idx,2)))*100));
            end
            if k == size(idx,2)
                fprintf('\n');
            end
        end
            
        end
        
        for i = 1:size(freq_mat1,1)
            for j = 1:size(freq_mat1,2)
                
                freq_mat1cell{i,j} = [freq_mat1cell{i,j} freq_mat1(i,j)];
                freq_mat2cell{i,j} = [freq_mat2cell{i,j} freq_mat2(i,j)];
                
            end
        end
    end
    
    freq_mat1cell = freq_mat1cell(1:size(freq_mat1,1), 1:size(freq_mat1,2));
    freq_mat2cell = freq_mat2cell(1:size(freq_mat1,1), 1:size(freq_mat1,2));
    
    
    freq_mat1 = cellfun(@mean, freq_mat1cell);
    freq_mat2 = cellfun(@mean, freq_mat2cell);
    
    maxfreq_store1 = mean(maxfreq_store1,1);
    maxfreq_store2 = mean(maxfreq_store2,1);
    %maxfreq_store3 = mean(maxfreq_store3,1);
    
    
end


% Plot

% 4 plot options: 
% 1) Heat plot spectrum, two separate plots
% 2) Heat plot difference in two spectrums, single plot
% 3) Line plot maxfreq varaition, two separate plots
% 4) Line plot difference in maxfreq, single plot

if plt ~= 0

time_vec = zeros(size(idx,2),1);
time_vec(1) = win_length/2;
    for i = 2:length(time_vec)
                  
            time_vec(i) = time_vec(i-1)+(win_length-overlap);
    end

if (plt == 1) || (plt == 2)
    
    plot_peak_freq_tf(freq_mat1, freq_mat2,time_vec, plt,...
        struct1_name, struct2_name, dataname, f1);
else
    plot_peak_freq_tf(maxfreq_store1, maxfreq_store2,time_vec,plt,...
        struct1_name, struct2_name, dataname, f1);
end

% TO USE THIS MUST CHANGE plot_peak_freq_tf.m ALSO!!!
% if (plt == 4)
%     plot_peak_freq_tf(maxfreq_store3, maxfreq_store2,time_vec,plt,...
%         struct1_name, struct2_name, dataname, f1);
% end
end
