function result = AnalyzeData(src,funcIn,varargin)
%ANALYZEDATA - Apply an analysis function to a DynaSim study or DS data, optionally saving data
%
% Pass a single DynaSim data structure or an array of data structures to a
% user-specified analysis function, add varied info to the results and
% optionally save the output structure.
%
% Usage:
%   result = AnalyzeData(data,func,'option1',value1,...) % pass data or datafile name
%   result = AnalyzeData(studyinfo,func,'option1',value1,...) % pass studyinfo struct
%   result = AnalyzeData(study_dir,func,'option1',value1,...) % pass study_dir containing studyinfo.mat
%
% Inputs:
%   - First inputs/argument:
%     - data: DynaSim data structure or data file path (s)
%     - studyinfo: DynaSim studyinfo structure or path to studyinfo
%     - study_dir: DynaSim study directory containing studyinfo.mat
%   - func: function handle or cell array of function handles pointing to plot
%           or analysis function(s)
%   - options: (key/value pairs are passed on to the analysis function)
%     'save_results_flag'   : whether to save result {0 or 1} (default: 0)
%     'result_file'         : where to save result (default: 'result.mat')
%     'format'              : format for saved plots if figures are generated
%                             {'svg','jpg','eps','png'} (default: 'svg')
%     'varied_filename_flag': whether to make filename based on the varied
%                             parameters and type of plot {0 or 1} (default: 0)
%     'function_options'    : cell array of option cell arrays {'option1',value1,...}
%                             in which each cell corresponds to the options for
%                             the corresponding function cell. if only passing a
%                             single func, can specificy function options as
%                             key,val list as varargin for AnalyzeData
%     'one_subplot_flag'    : whether to plot each set of data in a separate
%                             file. Make sure the plot functions include 'plot' 
%                             in the function name {0 or 1} (default: 1)
%
% Outputs:
%   - result: structure returned by the analysis function
%
% TODO: annotate figures with data set-specific modifications
%
% TODO: handle multiple plottypes as cell array, if postSimBool don't load all
%       data at once
%
% See also: SimulateModel

% check inputs
options = CheckOptions(varargin,{...
  'result_file','result',[],...
  'save_results_flag',0,{0,1},...
  'format','svg',{'svg','jpg','eps','png','fig'},...
  'varied_filename_flag',0,{0,1},...
  'plot_type','waveform',{'waveform','rastergram','raster','power','rates'},...
  'save_prefix',[],[],...
  'function_options',{},[],...
  'one_subplot_flag',1,{0,1},...
  },false);

% save data if no output is requested
if nargout<1
  options.save_results_flag = 1;
end

% determine type of src
if ischar(src)
  if exist(src,'file') % data file or studyinfo.mat
    if strfind(src, 'studyinfo') %studyinfo.mat
      [data,studyinfo] = ImportData(src, varargin{:}); % load data
      studyinfo.study_dir = fileparts(src);
    else % data file
      [data,studyinfo] = ImportData(src, varargin{:}); % load data
    end
  elseif isdir(src) % study_dir
    [data,studyinfo] = ImportData(src, varargin{:}); % load data
    studyinfo.study_dir = src;
  else
    try
      [data,studyinfo] = ImportData(src, varargin{:}); % load data
    catch
      error('Unknown source for first input/argument.')
    end
  end
elseif isstruct(src) && length(src)>1 % data file cell array
  data = src;
elseif isstruct(src) % single data struct or studyinfo struct
  if isfield(src,'time') % single data file 
    data = src;
  else % studyinfo struct
    [data,studyinfo] = ImportData(src, varargin{:}); % load data
  end
elseif iscell(src) % cell array of files
  [data,studyinfo] = ImportData(src, varargin{:}); % load data
else
  try
    [data,studyinfo] = ImportData(src, varargin{:}); % load data
  catch
    error('Unknown source for first input/argument.')
  end
end

% check study_dir
if ~isfield(studyinfo,'study_dir') || isempty(studyinfo.study_dir) || ~isdir(studyinfo.study_dir)
  studyinfo.study_dir = pwd;
end

% convert data to double precision before analysis
for j = 1:length(data)
  for k = 1:length(data(j).labels)
    fld = data(j).labels{k};
    data(j).(fld) = double(data(j).(fld));
  end
end

% check if func not a cell array
if isa(funcIn,'function_handle')
  nFunc = 1;
  funcIn = {funcIn}; % make into cell array
elseif ischar(funcIn)
  funcIn = {str2func(funcIn)}; % convert string to fn handle
  
  if ~isfunction(funcIn)
    error('Post-processing function must be supplied as a function handle or function name string');
  end
  
  nFunc = 1;
else
  nFunc = numel(funcIn);
end

% check if postSim
postSimBool = (length(data) > 1);

for fInd = 1:nFunc % loop over function inputs
  func = funcIn{fInd};
  
  % confirm function handle or convert to one if possible
  if ~isa(func,'function_handle')
    if ischar(func)
      func = str2func(func); % convert string to fn handle
      
      if ~isfunction(func)
        error('Post-processing function must be supplied as a function handle or function name string');
      end
    else
      error('Post-processing function must be supplied as a function handle or function name string');
    end
  end
  
  % change result_file if varied_filename_flag
  if options.varied_filename_flag && isfield(data, 'varied')
    if isfield(options, 'save_prefix') && ~isempty(options.save_prefix)
      prefix = options.save_prefix;
    else
      if regexpi(func2str(func), 'plot')
        prefix = options.plot_type;
      else
        prefix = func2str(func);
      end
    end

    options.result_file = nameFromVaried(data, prefix, options.result_file);
  end
  
  % do analysis
  fprintf('\tExecuting post-processing function: %s\n',func2str(func));
  tstart = tic;
  
  if ~(regexpi(func2str(func), 'plot') && options.one_subplot_flag) % skip if plot function and one_subplot_flag
    if isempty(options.function_options)
      result = feval(func,data,varargin{:});
    else
      function_options = options.function_options{fInd};
      result = feval(func,data,function_options{:});
    end
  end
  
  fprintf('\t\tElapsed time: %g sec\n',toc(tstart));
  
  % determine if result is a plot handle or derived data
  if options.one_subplot_flag || all(ishandle(result)) % analysis function returned a graphics handle
    if ~options.one_subplot_flag
      nResults = length(result);
    else
      nResults = length(data);
    end
    
    for iResult = 1:nResults
      % save plot
      if options.save_results_flag
        extension = ['.' options.format]; % '.svg'; % {.jpg,.svg}
        
        if nResults == 1
          fname = [options.result_file extension];
        elseif ~isempty(studyinfo)
          simID = studyinfo.simulations(iResult).sim_id;
          prefix = func2str(func);
          fname = [prefix '_sim' num2str(simID) '_plot' num2str(fInd) '_' func2str(func)];
          
          if options.one_subplot_flag
            result = feval(func,data(iResult),options.function_options{:});
            iResult = 1; % HACK
          end
          
          % change result_file if varied_filename_flag
          if options.varied_filename_flag && isfield(data, 'varied')
            if isfield(options, 'save_prefix') && ~isempty(options.save_prefix)
              prefix = options.save_prefix;
            else
              if regexpi(func2str(func), 'plot')
                plot_options = CheckOptions(varargin,{'plot_type',['plot' num2str(fInd)],{'waveform','rastergram','raster','power','rates'}},false);
                prefix = plot_options.plot_type;
              else
                prefix = func2str(func);
              end
            end

            fname = nameFromVaried(data, prefix, fname);
            
            fDir = fullfile(studyinfo.study_dir, 'postSimPlots');
            if ~exist(fDir,'dir')
              mkdir(fDir)
            end
            fPath = fullfile(fDir,fname);
            fprintf('\t\tSaving derived data: %s\n', fPath);
          end % varied_filename_flag
        else
          fname = [options.result_file '_page' num2str(iResult) extension];
        end
        
        set(gcf,'PaperPositionMode','auto');
        fprintf('\t\tSaving plot: %s\n',fname);
        
        switch extension
          case '.svg'
            plot2svg(fname,result(iResult));
          case '.jpg'
            print(result(iResult),fname,'-djpeg');
          case '.eps'
            print(result(iResult),fname,'-depsc');
          case '.png'
            print(result(iResult),fname,'-dpng');
          case '.fig'
            savefig(result(iResult),fname);
        end
        
        if postSimBool
          close(result)
        end
      end
    end
  else % analysis function returned derived data
    if isstruct(result)
      result = add_modifications(result);
      
      for iResult = 1:length(result)
        % add options to result structure
        if length(varargin)>1
          for j = 1:2:length(varargin)
            result(iResult).options.(varargin{j}) = varargin{j+1};
          end
        else
          result(iResult).options = [];
        end
      end %iResult
    end %isstruct
    
    % save derived data
    if options.save_results_flag
      if length(result) == 1
        fname = options.result_file;
        extension = '.mat';
        
        if ~strcmp(fname(end-3:end), extension) %check for .mat extension
          fname = [fname extension];
        end
        
        fprintf('\t\tSaving derived data: %s\n', fname);
        save(fname,'result','-v7.3');
      elseif ~isempty(studyinfo) % analyze study
        allResults = result;
        clear result;
        for iResult = 1:length(allResults)
          result = allResults(iResult);
          simID = studyinfo.simulations(iResult).sim_id;
          prefix = func2str(func);
          fname = [prefix '_sim' num2str(simID) '_analysis' num2str(fInd) '_' func2str(func) '.mat'];
          
          % change result_file if varied_filename_flag
          if options.varied_filename_flag && isfield(data, 'varied')
            if isfield(options, 'save_prefix') && ~isempty(options.save_prefix)
              prefix = options.save_prefix;
            else
                prefix = func2str(func);
            end

            fname = nameFromVaried(data, prefix, fname);
          end % varied_filename_flag
          
          fDir = fullfile(studyinfo.study_dir, 'analyzedData');
          if ~exist(fDir,'dir')
            mkdir(fDir)
          end
          fPath = fullfile(fDir,fname);
          fprintf('\t\tSaving derived data: %s\n', fPath);
          save(fname,'result','-v7.3');
        end %iResult
      end % if length(result)==1
    end % save_results_flag
  end % ishandle(result)
end % fInd

%% Nested functions
  function result = add_modifications(result)
    % add modifications to result structure, excluding modifications made
    % within experiments. note: while this nested function is similar to 
    % prepare_varied_metadata in SimulateModel, the data structure contains
    % all modifications (those within and across experiments; listed in 'varied'). 
    % the result structure collapses data sets from an experiment into a single
    % result; thus, each result corresponds to modifications across
    % experiments but not within them; those modifications are stored in
    % the simulator options.
    if ~isempty(data(1).simulator_options.modifications)
      varied = {};
      mods = data(1).simulator_options.modifications;
      for ii = 1:length(result)
        for jj = 1:size(mods,1) 
          % prepare valid field name for thing varied:
          fld = [mods{jj,1} '_' mods{jj,2}];
          
          % convert arrows and periods to underscores
          fld = regexprep(fld,'(->)|(<-)|(-)|(\.)','_');
          
          % remove brackets and parentheses
          fld = regexprep(fld,'[\[\]\(\)\{\}]','');
          result(ii).(fld) = mods{jj,3};
          varied{end+1} = fld;
        end
        result(ii).varied = varied;
        result(ii).modifications = mods;
      end
    elseif isfield(data,'varied') && length(data) == 1
      % add 'varied' info from data to result structure
      for ii = 1:length(result)
        result(ii).varied = data(1).varied;
        for jj = 1:length(data(1).varied)
          result(ii).(data(1).varied{jj}) = data(1).(data(1).varied{jj});
        end
      end
    end
  end % add_modifications

end %main fn
