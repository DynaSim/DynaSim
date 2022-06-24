function plot_peak_freq_tf(freqmat1, freqmat2, time_vec, plt,...
    structure1, structure2, dataset, freq_vec)
%function plot_peak_freq_tf(freqmat1, freqmat2, time_vec, plt,...
%    structure1, structure2, dataset, freq_vec)
%
% This function plots the time course of the power spectrum of two signals
% in a variety of different ways:
% 
% 1) Heat plot spectrum, two separate plots
% 2) Heat plot difference in two spectrums, single plot
% 3) Line plot maxfreq variation, two separate plots
% 4) Line plot difference in maxfreq, single plot
%
%
% INPUTS:
%
% freqmat1 - a matrix of either power spectrum values or a vector of 
% peak frequency values calculated over a series of time windows for
% sig_struct1
% freqmat2 - a matrix of either power spectrum values or a vector of 
% peak frequency values calculated over a series of time windows for
% sig_struct2
% time_vec - vector of centre times for windows along x-axis
% plt - either 1,2,3 or 4 depending on the plot option desired
% structure1 - (string) name of the brain structure1
% structure2 - (string) name of the brain structure2
% dataset - (string) name of the dataset used 
% freq_vec - vector of frequencies to plot on y-axis
%
% Author: Angela Onslow, 2011


figure

if plt == 1
    
    subplot(2,1,1)
    imagesc(time_vec,freq_vec, freqmat1)
    axis xy
    title(['Spectrogram:',' Dataset = ', dataset, ' Structure = ',...
        structure1])
    ylabel('Frequency / Hz')
    xlabel('Time / s ')
    colorbar
    
    subplot(2,1,2)
    imagesc(time_vec,freq_vec, freqmat2)
    axis xy
    title(['Spectrogram:',' Dataset = ', dataset, ' Structure = ',...
        structure2])
    ylabel('Frequency / Hz')
    xlabel('Time / s ')
    colorbar
end

if plt == 2
    imagesc(time_vec,freq_vec, (freqmat1-freqmat2))
    axis xy
    title(['Difference spectrogram:', 'Dataset = ', dataset, ' Structures = ',...
        structure1,'-',structure2])
    ylabel('Frequency / Hz')
    xlabel('Time / s ')
    colorbar
end

if plt == 3
    subplot(2,1,1)
    plot(time_vec,freqmat1)
    axis xy
    title(['Peak freq variation:',' Dataset = ', dataset, ' Structure = ',...
        structure1])
    ylabel('Frequency / Hz')
    xlabel('Time / s ')
    
    
    subplot(2,1,2)
    plot(time_vec,freqmat2)
    title(['Peak freq variation:',' Dataset = ', dataset, ' Structure = ',...
        structure2])
    ylabel('Frequency / Hz')
    xlabel('Time / s ')
end
  
if plt == 4
 
    plot(time_vec,(freqmat1-freqmat2))
    title(['Peak freq variation:', 'Dataset = ', dataset, ' Structures = ',...
        structure1,'-',structure2])
    ylabel('Frequency / Hz')
    xlabel('Time / s ')
end