function pac_plot_fun_tf(pacmat, time_vec, amp_freq_vec, measure,...
    structure1, structure2, dataset, ph_freq)
% function pac_plot_fun_tf(pacmat, time_vec, amp_freq_vec, measure,...
%   structure1, structure2, dataset, ph_freq)
%
% This function plots a PACgram with time along the x-axis and modulated 
% (higher) frequency along the y-axis
%
% INPUTS:
%
% pacmat - matrix of PAC values to be plotted
% time_vec - vector of centre times for windows along x-axis
% amp_freq_vec - vector of centre frequency values which the higher,
% modulated signal has been filtered at
% measure - PAC measure which has been calculated to generate pacmat,
% either 'esc', 'mi' or 'cfc'
% structure1 - (string) name of the brain structure recorded from to obtain
% the PAC containing signal 'sig1', displayed on y-axis
% structure2 - (string) name of the brain structure recorded from to obtain
% the modulating signal 'sig2', displayed on x-axis
% dataset - (string) name of the dataset used to reference 'sig1' and
%
% Author: Angela Onslow, 2011

figure

imagesc(time_vec, amp_freq_vec, pacmat)


axis xy
title(['PAC measure = ',measure,' Dataset = ', dataset, ' Structures = ',...
    structure2,'-',structure1, 'Modulating freq = ', num2str(ph_freq)])
ylabel(['Amp freq / Hz - ', structure1])
xlabel('Time / s ')
colorbar
