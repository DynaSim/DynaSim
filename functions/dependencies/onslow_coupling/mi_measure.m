 function [mival, mi_angle, mi_comodulogram, modulation_index] = mi_measure(phase_sig, amp_sig, calc_comodulograms, Fs)
% function mival = mi_measure(phase_sig, amp_sig)
%
% Returns a value for the MI measure calculated between two signals.
% (Functionality to deal with multiple trials will be added soon)
%
% INPUTS:
%
% phase_sig - the instantaneous phase values for a signal which has been
% filtered for a lower, modulating frequency, passed as a column vector
%
% amp_sig - the amplitude values for a signal which has been filtered for a
% higher, modulated frequency, passed as a column vector 
%
% Author: Angela Onslow, May 2010
% Edited David Stanley 2015 at Boston University, to work with multiple trials

plot_on = 0;

num_trials = size(phase_sig, 2);

for count = 1:num_trials
    
    %Create composite signal
    z = amp_sig(:,count).*exp(1i*phase_sig(:,count));
    m_raw(count) = mean(z);  %Compute the mean length of composite signal.
        

    mival(count,1) = ((m_raw(count)));
    
    
    if plot_on
        %%
        figure; plot(real(z),imag(z),'.');
        hold on; plot(real(m_raw(count)),imag(m_raw(count)),'rx');
    end
    
    if calc_comodulograms
        % AES
        parms.window_length = 2.0;
%         parms.window_length = 3.0;
        
        parms.window_overlap = 1.0;
        parms.number_bins = 12;

        % AES 20221127 MAMA MIA
        % Fs = 1000/0.01 / 10; % don't forget downsampling

        %% Get angle and amplitude from the respective signals via Hilbert Transform
        % Kramer's GLMCFC code is espcially concise/informative. Note the '.Data'
        %     phi = angle(hilbert(slow_data.Data));
        %     amp = abs(hilbert(fast_data.Data));
        phi = phase_sig;
        amp = amp_sig;

        %% Construct windows across the time series
        %    Adapted & taken from Angela Onslow's 'pac_code_best/window_data.m' of
        %    her MATLAB Toolbox for Estimating Phase-Amplitude Coupling from
        %
        %    http://www.cs.bris.ac.uk/Research/MachineLearning/pac/
        %   
        % Compute indices
        number_amp = size(amp,1);
        number_trials = size(amp,2);
        number_windows = ceil(parms.window_length*Fs);
        number_overlaps = ceil(parms.window_overlap*Fs);
        idx = bsxfun(@plus, (1:number_windows)', 1+(0:(fix((number_amp-number_overlaps)/(number_windows-number_overlaps))-1))*(number_windows-number_overlaps))-1;
    
        % Initialize the main data objects
        modulation_index_timeseries = [];
        mi_comodulogram = [];
        %% Loop over sliding windows
        for k=1:size(idx,2)
            amp_window = [];
            phi_window = [];
            % Loop over trials
            for j = 1:number_trials
                amp_window = [amp_window, amp(idx(:,k),j)];
                phi_window = [phi_window, phi(idx(:,k),j)];
            end
            %% Bin the faster frequency's amplitude in the slower's phase bins
            %    Adapted & taken from Richardson Leao's copy of Adriano Tort's
            %    'Neurodynamics-master/16ch/Comodulation/ModIndex_v1.m' of the
            %    'Neurodynamics-Toolbox' repo on Github, at
            %
            %    https://github.com/cineguerrilha/Neurodynamics
            %
            phi_bin_beginnings = zeros(1,parms.number_bins); % this variable will get the beginning (not the center) of each bin (in rads)
            bin_size = 2*pi/parms.number_bins;
            for j=1:parms.number_bins
                phi_bin_beginnings(j) = -pi+(j-1)*bin_size;
            end

            % Now we compute the mean amplicomodtude in each phase:
            amp_means = zeros(1,parms.number_bins);
            for j=1:parms.number_bins
                phi_indices = find((phi_window >= phi_bin_beginnings(j)) & (phi_window < phi_bin_beginnings(j)+bin_size));
                amp_means(j) = mean(amp_window(phi_indices));
            end
            mi_comodulogram = [mi_comodulogram, (amp_means/sum(amp_means))'];
    
            % Quantify the amount of amp modulation by means of a normalized entropy index (Tort et al PNAS 2008):
            modulation_index=(log(parms.number_bins)-(-sum((amp_means/sum(amp_means)).*log((amp_means/sum(amp_means))))))/log(parms.number_bins);
            modulation_index_timeseries = [modulation_index_timeseries, modulation_index];

%         % Debug, for mid-function plotting:
%         % So note that the center of each bin (for plotting purposes) is phi_bin_beginnings+bin_size/2
%         % at this point you might want to plot the result to see if there's any amplitude modulation
%         figure(10)
%         bar(10:10:360,(amp_means/sum(amp_means)),'k')
%         xlim([0 360])
%         set(gca,'xtick',0:180:360)
%         xlabel('Phase (Deg)')
%         ylabel('Amplitude')
        end
    else
        mi_comodulogram = 0;
        %% Bin the faster frequency's amplitude in the slower's phase bins
        %    Adapted & taken from Richardson Leao's copy of Adriano Tort's
        %    'Neurodynamics-master/16ch/Comodulation/ModIndex_v1.m' of the
        %    'Neurodynamics-Toolbox' repo on Github, at
        %
        %    https://github.com/cineguerrilha/Neurodynamics
        %
        nbins = 18;
        phase_bin_beginnings = zeros(1,nbins); % this variable will get the beginning (not the center) of each bin (in rads)
        bin_size = 2*pi/nbins;
        for j=1:nbins
            phase_bin_beginnings(j) = -pi+(j-1)*bin_size;
        end

        % Now we compute the mean amplitude in each phase:
        amp_means = zeros(1,nbins);
        for j=1:nbins
            phase_indices = find((phase_sig >= phase_bin_beginnings(j)) & (phase_sig < phase_bin_beginnings(j)+bin_size));
            amp_means(j) = mean(amp_sig(phase_indices));
        end
        
        mi_comodulogram = (amp_means/sum(amp_means))';
    
        % Quantify the amount of amp modulation by means of a normalized entropy index (Tort et al PNAS 2008):
        % Note: If the MI is NaN, then likely one of the amp_means is a NaN
        % owing to one of the 18 phase bins not being detected at all; this
        % clearly indicates that our SWO-filtered phase signal is missing
        % at least 1/18th of the cycle, and therefore this
        % particularly-filtered signal should probably be thrown out.
        modulation_index=(log(nbins)-(-sum((amp_means/sum(amp_means)).*log((amp_means/sum(amp_means))))))/log(nbins);
    end
end

if num_trials > 1
    mival = mean(mival);
end

mival = abs(mival);     % Dave's edit - moved the absolute value to outside
                        % the loop to work with multiple trials.
                        % This requires that the phase of coupling to be
                        % consistent across trials. (i.e. the gamma
                        % oscillations must always appear at the same phase
                        % of the theta peaks.)

mi_angle = angle(m_raw);

end
