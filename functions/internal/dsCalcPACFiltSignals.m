function data = dsCalcPACFiltSignals(data, varargin)
%CALCPOWER - Compute spectral analysis of DynaSim data
%
% Usage:
%   data = dsCalcPower(data,'option',value)
%
% Inputs:
%   - data: DynaSim data structure (see dsCheckData)
%   - options:
%     'variable'              : name of field containing data on which to
%                               calculate firing rates (default: *_spikes or
%                               first variable in data.labels)
%     'time_limits'           : [beg,end] (units of data.time)
%     'smooth_factor'         : number of samples for smoothing the spectrum (default: 5)
%                               - tip: set to 1 to avoid smoothing
%   - options for peak detection:
%     'min_peak_frequency'    : Hz, min frequency for peak detection (default: 2)
%     'max_peak_frequency'    : Hz, max frequency for peak detection (default: 150)
%     'peak_threshold_prctile': percentile for setting power threshold for peak
%                               detection (default: 95)
%     'peak_area_width'       : Hz, size of frequency bin (centered on peak)
%                               over which to calculate area under spectrum (default: 5)
%     'exclude_data_flag'     : whether to remove simulated data from result
%                               structure (default: 0)
%
% Outputs:
%   - data structure organization:
%     data.VARIABLE_Power_SUA.frequency: TODO
%     data.VARIABLE_Power_SUA.PeakArea: area under spectrum around peak (one value per cell)
%     data.VARIABLE_Power_SUA.PeakFreq: frequency of spectral power (one value per cell)
%     data.VARIABLE_Power_SUA.Pxx: spectral power
%   - if populations present, data also includes:
%     data.VARIABLE_Power_MUA.frequency: TODO
%     data.VARIABLE_Power_MUA.PeakArea: TODO
%     data.VARIABLE_Power_MUA.PeakFreq: TODO
%     data.VARIABLE_Power_MUA.Pxx: spectrum of the mean waveform
%       - population mean spectrum of the individual waveforms can be
%           calculated using "mean(data.VARIABLE_Power_MUA.Pxx,2)".
%   - Note:
%     - "VARIABLE" can be specified as the name of a variable listed in
%         data.labels, a cell array of string listing variable names, or as a
%         regular expression pattern for identifying variables to process. See
%         dsSelectVariables for more info on supported specifications.
%
% Examples:
%   s=[];
%   s.populations(1).name='E';
%   s.populations(1).equations='dv[2]/dt=@current+10; {iNa,iK}; v(0)=-65';
%   s.populations(2).name='I';
%   s.populations(2).equations='dv/dt=@current+10; {iNa,iK}; v(0)=-65';
%   data=dsSimulate(s,'tspan',[0 1000]);
%   data=dsCalcPower(data,'variable','v');
%   % Plot the spectrum of the E-cell average population voltage
%   figure; plot(data.E_v_Power_MUA.frequency,data.E_v_Power_MUA.Pxx);
%   xlabel('frequency (Hz)'); ylabel('power'); xlim([0 200]);
%
% See also: PlotPower, dsAnalyzeStudy, dsSimulate, dsCheckData, dsSelectVariables
%
% Author: Jason Sherfey, PhD <jssherfey@gmail.com>
% Copyright (C) 2016 Jason Sherfey, Boston University, USA

%% 1.0 Check inputs
options=dsCheckOptions(varargin,{...
  'variable',[],[],...
  'time_limits',[-inf inf],[],...
  'smooth_factor',5,[],... % number of samples for smoothing the spectrum
  'min_peak_frequency',1,[],... % Hz, min frequency for peak detection
  'max_peak_frequency',200,[],... % Hz, max frequency for peak detection
  'peak_threshold_prctile',95,[],... % percentile for setting power threshold for peak detection
  'peak_area_width',5,[],... % Hz, size of frequency bin (centered on peak) over which to calculate area under spectrum
  'exclude_data_flag',0,{0,1},...
  'timeBandwidthProduct',[],[],... % time-bandwidth product for multi-taper method
  'output_suffix','',[],...
  'auto_gen_test_data_flag',0,{0,1},...
  'phase_freqs', [0.5:0.5:2.0],[],...% Hz, Frequencies to analyze for the phase of the "slow" or modulating signal for comodulograms
  'ampl_freqs', [8:3:14],[],...% Hz, Frequencies to analyze for the amplitude of the "fast" or carrier signal for comodulograms
  'measure', 'mi',{'mi','esc','cfc'},...% Type of coupling measure
  'plt', 'n',[],...% Don't use internal code to plot data
  'waitbar', 0,[],...% Don't print ongoing progress of significance analysis
  'width', 7,[],...% Width of Morlet wavelets to use for filtering, whatever?
  'nfft', 2500000,[],... % AES TODO Samples to use for each time/freq bin
  'num_shf', 0,[],...% Don't run any statistical significance analysis on coupling
  'calc_comodulograms', 1, [],...
},false);


%% 'nfft', ceil(Fs/(diff(phase_freqs(1:2)))),[],...

data = dsCheckData(data, varargin{:});
% note: calling dsCheckData() at beginning enables analysis function to
% accept data matrix [time x cells] in addition to DynaSim data structure.

if numel(data)>1
  % use dsAnalyzeStudy to recursively call dsCalcPower on each data set
  data=dsAnalyzeStudy(data,@dsCalcPACFiltSignals,varargin{:});
  return;
end

% time parameters
time = data.time; % time vector
dt = time(2)-time(1); % time step
%% % ntime=length(time); % number of time points in full data set
t1=nearest(time,options.time_limits(1)); % index to first sample
t2=nearest(time,options.time_limits(2)); % index to last sample
nsamp=t2-t1+1; % number of samples for spectral estimate

% frequency parameters
Fs = fix(1/(dt/1000)); % effective sampling frequency

%% 2.0 set list of variables to process as cell array of strings
options.variable=dsSelectVariables(data(1).labels,options.variable, varargin{:});

%% 3.0 calculate power spectrum for each variable
if ~isfield(data,'results')
  data.results={};
end

warning off
for v=1:length(options.variable)
  % extract this data set
  var=options.variable{v};
  dat=data.(var);

  % determine how many cells are in this data set
  ncells=size(dat,2);

    % calculate MUA
    % TODO AES try with and without detrend
    X=detrend(nanmean(dat(t1:t2,:),2)); % detrend the data

    %% [tmpPxx,f] = pwelch(X,NFFT,[],NFFT,Fs); % calculate power

    fprintf('About to start running phase-amplitude coupling / comodulogram analysis\n')
    [filtered_slow, filtered_fast] = ...
        find_pac_shf_save_sigs(X, Fs, options.measure, X, ...
            options.phase_freqs, options.ampl_freqs, options.plt,...
            options.waitbar, options.width, options.nfft, ...
            options.num_shf, options.calc_comodulograms);

    %% Add resulting power spectra to data structure
    % organization scheme:
    % data.VARIABLE_Power_SUA.(Pxx,PeakFreq,PeakArea,frequency)
    % data.VARIABLE_Power_MUA.(Pxx,PeakFreq,PeakArea,frequency)
    data.([var '_FiltSigs' options.output_suffix]).filtered_slow=filtered_slow;
    data.([var '_FiltSigs' options.output_suffix]).filtered_fast=filtered_fast;
%     % AES debug
%     figure(10)
%     for ii = 1:size(comodulograms, 1)
%         for jj = 1:size(comodulograms, 2)
%             subplot(size(comodulograms,1), size(comodulograms,2), ii*jj)
%             imagesc(comodulograms{ii,jj})
%         end
%     end
%     
    if ~ismember([var '_FiltSign' options.output_suffix],data.results)
      data.results{end+1}=[var '_FiltSigs' options.output_suffix];
    end

    if options.exclude_data_flag
      for l=1:length(data.labels)
        data=rmfield(data,data.labels{l});
      end
    end
end % end of options.variable loop
