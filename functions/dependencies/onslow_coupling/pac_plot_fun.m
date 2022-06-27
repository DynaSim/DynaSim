function pac_plot_fun(pacmat, ph_freq_vec, amp_freq_vec, measure,...
    structure1, structure2, dataset, pac_angles)
% function pac_plot_fun(pacmat, ph_freq_vec, amp_freq_vec, measure,...
%   structure1, structure2, dataset)
%
% This function plots a PACgram.
%
% INPUTS:
%
% pacmat - matrix of PAC values to be plotted
% ph_freq_vec - vector of centre frequency values which the lower,
% modulating signal has been filtered at
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
% Author: Angela Onslow, 2010

figure

if strcmp(measure, 'esc')
    imagesc(ph_freq_vec, amp_freq_vec, pacmat)
else
    imagesc(ph_freq_vec, amp_freq_vec, pacmat)
end

hold on
quiver(ph_freq_vec, amp_freq_vec, ones(size(pac_angles)), pac_angles,...
    'Color', 'k', 'LineWidth', 2.0, 'AutoScaleFactor', 0.5)
hold off

axis xy
title(['PAC measure = ',measure,' Dataset = ', dataset, ' Structures = ',...
    structure2,'-',structure1])
ylabel(['Amp freq / Hz - ', structure1])
xlabel(['Phase freq / Hz - ', structure2])
colorbar




