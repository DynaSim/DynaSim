function pac_plot_fun_spec(pacmat, ph_freq_vec, amp_freq_vec, measure, ...
    structure1, structure2, dataset, sig1, Fs, sig2)
% function pac_plot_fun_spec(pacmat, ph_freq_vec, amp_freq_vec, measure, ...
%   structure1, structure2, dataset, sig1, Fs, sig2)
%
% This function plots a PACgram, with the frequency spectra of the two
% signals under comparison as a subplot below.
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
% 'sig2', displayed in title of plot
% sig1 - signal containing (or believed to contain) PAC
% Fs - sampling frequency, in Hz
% sig2 - signal containing (or believed to contain) the modulating signal
%
% Author: Angela Onslow, 2010



thesame = (sig1 == sig2);

figure

if thesame
    num_plots = 2;
else
    num_plots = 3;
end
subplot (num_plots, 1, 1);

if strcmp(measure, 'esc')
    imagesc(ph_freq_vec, amp_freq_vec, pacmat)
else
    imagesc(ph_freq_vec, amp_freq_vec, pacmat)
end
axis xy
title(['PAC measure = ',measure,' Dataset = ', dataset, ' Structures = ',...
    structure2,'-',structure1])
ylabel(['Amp freq / Hz - ', structure1])
xlabel(['Phase freq / Hz - ', structure2])
colorbar


[spec maxfreq specpeak f] =  freq_spec(mean(sig1, 2), Fs, 'n', max(amp_freq_vec));
%freq_spec(mean(sig1, 2), Fs, 'y', max(amp_freq_vec));
subplot(num_plots,1,2);
plot(f,spec)
title(['Frequency Spectrum -', structure1]);
xlabel('Frequency Hz')
ylabel('Amplitude')
h1 = subplot(num_plots,1,1);
h2 = subplot(num_plots,1,2);
p1 = get(h1, 'position');
p2 = get(h2, 'position');
set(h2,'Position',[p2(1) p2(2) p1(3) p1(4)]);

if ~thesame  
    [spec maxfreq specpeak f] =  freq_spec(mean(sig2, 2), Fs, 'n', max(amp_freq_vec));
    subplot(3,1,3)
    plot(f,spec)
    %freq_spec(mean(sig2, 2), Fs, 'y', max(amp_freq_vec));
    title(['Frequency Spectrum -', structure2]);
    xlabel('Frequency Hz')
    ylabel('Amplitude')
    h1 = subplot(3,1,1);
    h3 = subplot(3,1,3);
    p1 = get(h1, 'position');
    p3 = get(h3, 'position');
    set(h3,'Position',[p3(1) p3(2) p1(3) p1(4)]);
end


