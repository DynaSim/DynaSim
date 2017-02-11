function data = CalcPower(data,varargin)
%% data = CalcPower(data,'option',value)
% Inputs:
%   data - DynaSim data structure (see CheckData)
%   options:
%     'variable' - name of field containing data on which to calculate firing
%                rates (default: *_spikes or first variable in data.labels)
%     'time_limits' - [beg,end] (units of data.time)
%     'smooth_factor' - number of samples for smoothing the spectrum (default: 5)
%                       tip: set to 1 to avoid smoothing.
%   options for peak detection:
%     'min_peak_frequency' - Hz, min frequency for peak detection (default: 2)
%     'max_peak_frequency' - Hz, max frequency for peak detection (default: 150)
%     'peak_threshold_prctile' percentile for setting power threshold for peak detection (default: 95)
%     'peak_area_width' - Hz, size of frequency bin (centered on peak) over which to calculate area under spectrum (default: 5)
%     'exclude_data_flag' - whether to remove simulated data from result structure (default: 0)
% 
% Outputs:
%   data: data structure with spectral power in data.VARIABLE_Power_SUA.Pxx
%         data.VARIABLE_Power_SUA.PeakFreq: frequency of spectral power (one value per cell)
%         data.VARIABLE_Power_SUA.PeakArea: area under spectrum around peak (one value per cell)
%   NOTE: for populations: spectrum of the mean waveform is stored in 
%         data.VARIABLE_Power_MUA.Pxx. population mean spectrum of the individual
%         waveforms can be calculated as mean(data.VARIABLE_Power_MUA.Pxx,2).
% 
% organization scheme for spectral results:
% data.VARIABLE_Power_SUA.(Pxx,PeakFreq,PeakArea,frequency)
% data.VARIABLE_Power_MUA.(Pxx,PeakFreq,PeakArea,frequency)
% 
% note:
% "variable" can be specified as the name of a variable listed in
% data.labels, a cell array of string listing variable names, or as a 
% regular expression pattern for identifying variables to process.
% See SelectVariables for more info on supported specifications.
% 
% Examples:
% s=[];
% s.populations(1).name='E';
% s.populations(1).equations='dv[2]/dt=@current+10; {iNa,iK}; v(0)=-65';
% s.populations(2).name='I';
% s.populations(2).equations='dv/dt=@current+10; {iNa,iK}; v(0)=-65';
% data=SimulateModel(s,'tspan',[0 1000]);
% data=CalcPower(data,'variable','v');
% % Plot the spectrum of the E-cell average population voltage
% figure; plot(data.E_v_Power_MUA.frequency,data.E_v_Power_MUA.Pxx); 
% xlabel('frequency (Hz)'); ylabel('power'); xlim([0 200]);
% 
% See also: PlotPower, AnalyzeStudy, SimulateModel, CheckData, SelectVariables

%% 1.0 Check inputs
options=CheckOptions(varargin,{...
  'variable',[],[],...        
  'time_limits',[-inf inf],[],...
  'smooth_factor',5,[],... % number of samples for smoothing the spectrum
  'min_peak_frequency',2,[],... % Hz, min frequency for peak detection
  'max_peak_frequency',200,[],... % Hz, max frequency for peak detection
  'peak_threshold_prctile',95,[],... % percentile for setting power threshold for peak detection
  'peak_area_width',5,[],... % Hz, size of frequency bin (centered on peak) over which to calculate area under spectrum
  'exclude_data_flag',0,{0,1},...
  },false);

data = CheckData(data);
% note: calling CheckData() at beginning enables analysis function to
% accept data matrix [time x cells] in addition to DynaSim data structure.

if numel(data)>1
  % use AnalyzeStudy to recursively call CalcPower on each data set
  data=AnalyzeStudy(data,@CalcPower,varargin{:});
  return;
end

% time parameters
time = data.time; % time vector
dt = time(2)-time(1); % time step
ntime=length(time); % number of time points in full data set
t1=nearest(time,options.time_limits(1)); % index to first sample
t2=nearest(time,options.time_limits(2)); % index to last sample
nsamp=t2-t1+1; % number of samples for spectral estimate

% frequency parameters
Fs = fix(1/(dt/1000)); % effective sampling frequency
Fmin=options.min_peak_frequency; % min frequency for peak detection
Fmax=options.max_peak_frequency; % max frequency for peak detection
Fwin=options.peak_area_width; % size of frequency bin around peak for calculating area under spectrum
thresh_prctile=options.peak_threshold_prctile; % percentile for setting power threshold for peak detection
smooth_factor=options.smooth_factor; % number of samples to smooth spectrum
Fwin=options.peak_area_width; % size of frequency bin (centered on peak) over which to calculate area under spectrum
FreqRange=[max(Fmin,2/time(end)) Fmax]; % range to look for spectral peaks
NFFT=2^(nextpow2(nsamp-1)-1);%2); % <-- use higher resolution to capture STO freq variation
WINDOW=2^(nextpow2(NFFT-1)-3);
NOVERLAP=[]; % spectral parameters

%% 2.0 set list of variables to process as cell array of strings
options.variable=SelectVariables(data(1).labels,options.variable);

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
  % preallocation
  PeakFreq=nan(1,ncells);
  PeakArea=nan(1,ncells);   
  % SUA spectra: loop over cells
  for i=1:ncells
    % select data
    X=detrend(dat(t1:t2,i)); % detrend the data
    % calculate spectral estimate
    [tmpPxx,f] = pmtm(X, [], NFFT, Fs); % pwelch(X,NFFT,[],NFFT,Fs); % calculate power
    if i==1
      % get size of spectrum and preallocate result matrix
      nfreq=length(f);
      Pxx=nan(nfreq,ncells);
    end
    if all(isnan(tmpPxx(:)))
      tmpPxx=zeros(size(tmpPxx));
    end
    if ~isa(tmpPxx,'double')
      % convert to double precision
      tmpPxx=double(tmpPxx);
    end
    if smooth_factor>1
      % smooth the spectrum
      tmpPxx=smooth(tmpPxx,smooth_factor);
    end
    % Peak Detection:
    % select range of frequencies over which to look for peaks
    sel = find(FreqRange(1)<=f & f<=FreqRange(end));
    % set threshold for peak detection
    ht=prctile(log10(tmpPxx(sel)),thresh_prctile);
    if ~isnan(ht)
      % get index of peaks in range over threshold
      [PeakPower,PPind]=findpeaks(log10(tmpPxx(sel)),'MinPeakHeight',ht,'NPeaks',3);
    else
      PPind=[];
    end
    if ~isempty(PPind)
      % if multiple peaks, only consider the largest
      if numel(PPind)>1
        PPind=PPind(max(PeakPower)==PeakPower); %PPind=PPind(1);
      end
      % get frequency at that index
      PeakFreq(i) = f(sel(PPind));
      % set limits for calculating area under spectrum around peak
      flo=PeakFreq(i)-Fwin/2;
      fhi=PeakFreq(i)+Fwin/2;
      sel2=(flo<=f & f<=fhi);
      % calculate area under spectrum around peak
      PeakArea(i) = sum(tmpPxx(sel2))*(f(2)-f(1));
    else
      PeakFreq(i)=nan;
      PeakArea(i)=nan;
    end
    % Store results
    Pxx(:,i)=tmpPxx;
  end
  % -----------------------------------------------------
  % Repeat spectral estimate for MUA: 
  if ncells==1
    % same as SUA
    Pxx_mean=Pxx;
    Pxx_mean_PeakFreq=PeakFreq;
    Pxx_mean_PeakArea=PeakArea;
  else
    % calculate MUA
    X=detrend(nanmean(dat(t1:t2,:),2)); % detrend the data
    % calculate spectral estimate
    [tmpPxx,f] = pwelch(X,NFFT,[],NFFT,Fs); % calculate power
    if all(isnan(tmpPxx(:)))
      tmpPxx=zeros(size(tmpPxx));
    end
    if ~isa(tmpPxx,'double')
      % convert to double precision
      tmpPxx=double(tmpPxx);
    end
    if smooth_factor>1
      % smooth the spectrum
      tmpPxx=smooth(tmpPxx,smooth_factor);
    end
    % Peak Detection:
    % select range of frequencies over which to look for peaks
    sel = find(FreqRange(1)<=f & f<=FreqRange(end));
    % set threshold for peak detection
    ht=prctile(log10(tmpPxx(sel)),thresh_prctile);
    if ~isnan(ht)
      % get index of peaks in range over threshold
      [PeakPower,PPind]=findpeaks(log10(tmpPxx(sel)),'MinPeakHeight',ht,'NPeaks',3);
    else
      PPind=[];
    end    
    if ~isempty(PPind)
      % if multiple peaks, only consider the largest
      if numel(PPind)>1
        PPind=PPind(max(PeakPower)==PeakPower); %PPind=PPind(1);
      end
      % get frequency at that index
      Pxx_mean_PeakFreq = f(sel(PPind));
      % set limits for calculating area under spectrum around peak
      flo=Pxx_mean_PeakFreq-Fwin/2;
      fhi=Pxx_mean_PeakFreq+Fwin/2;
      sel2=(flo<=f & f<=fhi);
      % calculate area under spectrum around peak
      Pxx_mean_PeakArea = sum(tmpPxx(sel2))*(f(2)-f(1));
    else
      Pxx_mean_PeakFreq=nan;
      Pxx_mean_PeakArea=nan;
    end    
    Pxx_mean=tmpPxx;
  end
  
  % Add resulting power spectra to data structure
  % organization scheme:
  % data.VARIABLE_Power_SUA.(Pxx,PeakFreq,PeakArea,frequency)
  % data.VARIABLE_Power_MUA.(Pxx,PeakFreq,PeakArea,frequency)
  data.([var '_Power_SUA']).Pxx=Pxx;
  data.([var '_Power_SUA']).PeakFreq=PeakFreq;
  data.([var '_Power_SUA']).PeakArea=PeakArea;
  data.([var '_Power_SUA']).frequency=f;
  data.([var '_Power_MUA']).Pxx=Pxx_mean;
  data.([var '_Power_MUA']).PeakFreq=Pxx_mean_PeakFreq;
  data.([var '_Power_MUA']).PeakArea=Pxx_mean_PeakArea;
  data.([var '_Power_MUA']).frequency=f;
  if ~ismember([var '_Power_SUA'],data.results)
    data.results{end+1}=[var '_Power_SUA'];
  end
  if ~ismember([var '_Power_MUA'],data.results)
    data.results{end+1}=[var '_Power_MUA'];
  end
  if options.exclude_data_flag
    for l=1:length(data.labels)
      data=rmfield(data,data.labels{l});
    end
  end
%   % alternate organization scheme:
%   data.([var '_Pxx'])=Pxx;
%   data.([var '_Pxx_PeakFreq'])=PeakFreq;
%   data.([var '_Pxx_PeakArea'])=PeakArea;
%   data.([var '_Pxx_mean'])=Pxx_mean;
%   data.([var '_Pxx_mean_PeakFreq'])=Pxx_mean_PeakFreq;
%   data.([var '_Pxx_mean_PeakArea'])=Pxx_mean_PeakArea;
%   if ~ismember([var '_Pxx'],data.results)  
%     data.results{end+1}=[var '_Pxx'];
%     data.results{end+1}=[var '_Pxx_PeakFreq'];
%     data.results{end+1}=[var '_Pxx_PeakArea'];
%     data.results{end+1}=[var '_Pxx_mean'];
%     data.results{end+1}=[var '_Pxx_mean_Pxx_PeakFreq'];
%     data.results{end+1}=[var '_Pxx_mean_Pxx_PeakArea'];
%   end
  
end
