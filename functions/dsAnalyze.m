function result = dsAnalyze(data,func,varargin)
%DSANALYZE - Apply an analysis function to DynaSim data, optionally saving data
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
%           or analysis function(s). Should not contain case-insensitive string
%           'plot' unless is a function that returns a figure handle for plotting.
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
%     'load_all_data_flag'  : whether to load all the data in studyinfo
%                             at once {0 or 1} (default: 0)
%
% Outputs:
%   - result: structure returned by the analysis function
%
% TODO: annotate figures with data set-specific modifications
%
%
% See also: dsSimulate

%% General cases:
%   - data struct (likely from SimualteModel call)
%   - data struct array
%   - studyinfo with load_all_data_flag==0
%   - studyinfo with load_all_data_flag==1

% check inputs
options=ds.checkOptions(varargin,{...
  'result_file','result',[],...
  'save_results_flag',0,{0,1},...
  'format','svg',{'svg','jpg','eps','png','fig'},...
  'varied_filename_flag',0,{0,1},...
  'plot_type','waveform',{'waveform','rastergram','raster','power','rates','imagesc','heatmapFR','heatmap_sortedFR','meanFR','meanFRdens'},...
  'save_prefix',[],[],...
  'function_options',{},[],...
  'simIDs',[],[],...
  'load_all_data_flag',0,{0,1},...
  'xUnit',0,{0,1},...
  },false);

%% Setup unit testing
if options.xUnit
  fPaths = {};
  % shortcut for localfunction testing
  if isempty(src) && isempty(funcIn)
    result = struct('result',[], 'fPaths',fPaths, 'localfunctions',localfunctions);
    return
  end
  options.result_file = ds.nameFromVaried(data, prefix, options.result_file);
end


%% Save data if no output is requested.
if nargout<1
  options.save_results_flag = 1;
end

% varinputs = varargin;
% % check if simIDs specified
% if ~isempty(options.simIDs)
%   if isstruct(varinputs{1})
%     varinputs.simIDs = options.simIDs;
%   else
%     % find simIDs in varinputs
%     varinputs{find(~cellfun(@isempty,strfind(varinputs(1:2:end), 'simIDs')))*2} = options.simIDs;
%   end
% end

%% Parse src.
[data, studyinfo] = parseSrc(src, options, varargin);
% Data at this point:
%   - 'data' as single struct or struct array, or empty
%   - 'studyinfo' struct or empty

% check if study_dir defined
if isempty(studyinfo)
  studyinfoBool = false;
else
  studyinfoBool = true;
  if ~isfield(studyinfo,'study_dir') || isempty(studyinfo.study_dir) || ~isdir(studyinfo.study_dir)
    studyinfo.study_dir = pwd;
  end
end

% Data at this point:
%   - 'data' as single struct or struct array, or empty
%   - 'studyinfo' struct with many flds or just 'study_dir' field

% convert data to double precision before analysis
for j = 1:length(data)
  for k = 1:length(data(j).labels)
    fld = data(j).labels{k};
    data(j).(fld) = double(data(j).(fld));
  end
end

% convert func to handle if not a cell array
[funcIn, nFunc] = parseFuncIn(funcIn);

% check if postSim
postSimBool = studyinfoBool || (length(data) > 1); % since length(data)==1 and no studyinfo with SimulateModel call

for fInd = 1:nFunc % loop over function inputs
  func = funcIn{fInd};

  % confirm func is function handle or convert to one if possible
  func = parseFunc(func);

  % check if plot in fn name
  plotFnBool = ~isempty(regexpi(func2str(func), 'plot'));

  % change result_file if varied_filename_flag
  if options.varied_filename_flag && isfield(data, 'varied')
    options.result_file = filenameFromVaried(options.result_file, func, data, plotFnBool, options);
  end

  % do analysis
  fprintf('\tExecuting post-processing function: %s\n',func2str(func));
  tstart = tic;

  %% Eval func
  if length(data)==1 || postSimBool % Don't need to add check on load_all_data_flag, since if false data is empty.
    result = evalFnWithArgs(fInd, data, func, options, varargin);
  else
    result = [];
  end



  % calc nResults
  if ~isempty(results)
    nResults = length(result);
  elseif ~isempty(data)
    nResults = length(data);
  else
    nResults = length(siminfo.simulations);
  end

  fprintf('\t\tElapsed time: %g sec\n',toc(tstart));

  % Dave: Not all plotting functions will return a plot handle. For
  % example, dsPlot2 returns a nested structure of figure, axis, and plot
  % handles. This command updates it.
  if isstruct(result)
      if isfield(result,'hcurr')
          result = result.hcurr;
      end
  end

  % determine if result is a plot handle or derived data
  if all(ishandle(result)) || plotFnBool % analysis function returned a graphics handle or has plot in fn name
    %% Plot Function

    % will save plots else return main fn
    if options.save_results_flag
      % loop through results. all results may exist or need to be made during loop
      for iResult = 1:nResults
        extension = ['.' options.format]; % '.svg'; % {.jpg,.svg}

        if ~postSimBool % approx nResults == 1 && ~studyinfoBool % approx ~postSimBool
          fname = [options.result_file extension];
          fPath = fname;

          thisResult = result(iResult);
        elseif studyinfoBool
          simID = studyinfo.simulations(iResult).sim_id;
          prefix = func2str(func);
          fname = [prefix '_sim' num2str(simID) '_plot' num2str(fInd) '_' func2str(func)];

          if (postSimBool && ~options.load_all_data_flag)
            data = loadDataFromSingleSim(simID, options, varargin);

            %skip if no data
            if isempty(data)
              continue
            end

            % calc result for this data
            thisResult = evalFnWithArgs(fInd, data, func, options, varargin);
          end

          % change result_file if varied_filename_flag
          if options.varied_filename_flag && isfield(data, 'varied')
            fname = filenameFromVaried(fname, func, data, plotFnBool, options);
          end % varied_filename_flag

          % make fPath
          fDir = fullfile(studyinfo.study_dir, 'postSimPlots');
          if ~exist(fDir,'dir')
            mkdir(fDir)
          end
          fPath = fullfile(fDir,fname);
        else % length(result)>1 and ~studyinfoBool
          fname = [options.result_file '_page' num2str(iResult) extension];

          % make fPath
          fDir = fullfile(studyinfo.study_dir, 'postSimPlots');
          if ~exist(fDir,'dir')
            mkdir(fDir)
          end
          fPath = fullfile(fDir,fname);

          thisResult = result(iResult);
        end
        % Data needed for plotting:
        %   - thisResult
        %   - fPath

        set(thisResult, 'PaperPositionMode','auto');
        fprintf('\t\tSaving plot: %s\n',fname);

        if ~options.xUnit
          switch extension
            case '.svg'
              plot2svg(fpath,thisResult);
            case '.jpg'
              print(thisResult,fPath,'-djpeg');
            case '.eps'
              print(thisResult,fPath,'-depsc');
            case '.png'
              print(thisResult,fPath,'-dpng');
            case '.fig'
              savefig(thisResult,fPath);
          end
        else
          fPaths{end+1} = fPath;
        end

        if nResults > 1
          close(thisResult)
        end
      end %nResults
    end %save_results_flag
  else % analysis function returned derived data
    %% Analysis Function
    if isstruct(result)
      result = add_modifications(result, data);

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

    % save derived data else return main function
    if options.save_results_flag
      if studyinfoBool
        allResults = result;
        clear result;

        for iResult = 1:nResults
          simID = studyinfo.simulations(iResult).sim_id;

          if options.load_all_data_flag
            result = allResults(iResult);
          else % load data
            data = loadDataFromSingleSim(simID, options, varargin);

            %skip if no data
            if isempty(data)
              continue
            end

            % calc result for this data
            result = evalFnWithArgs(fInd, data, func, options, varargin);
          end

          prefix = func2str(func);
          fname = [prefix '_sim' num2str(simID) '_analysis' num2str(fInd) '_' func2str(func) '.mat'];

          % change result_file if varied_filename_flag
          if options.varied_filename_flag && isfield(data, 'varied')
            fname = filenameFromVaried(fname, func, data, plotFnBool, options);
          end % varied_filename_flag

          % make fPath
          fDir = fullfile(studyinfo.study_dir, 'analyzedData');
          if ~exist(fDir,'dir')
            mkdir(fDir)
          end
          fPath = fullfile(fDir,fname);

          fprintf('\t\tSaving derived data: %s\n', fPath);
          if ~options.xUnit
            save(fPath,'result','-v7.3');
          else
            fPaths{end+1} = fname;
          end
        end %iResult
      else % ~studyinfoBool, whether 1 or array of struct
        fname = options.result_file;
        extension = '.mat';

        if ~strcmp(fname(end-3:end), extension) %check for .mat extension
          fname = [fname extension];
        end

        fprintf('\t\tSaving derived data: %s\n', fname);
        if ~options.xUnit
          save(fname,'result','-v7.3');
        else
          fPaths{end+1} = fname;
        end
      end % scenarios
    end % save_results_flag
  end % ishandle(result)
end % fInd

%% Unit Testing
if options.xUnit
  result = struct('result',result, 'fPaths',fPaths, 'localfunctions',localfunctions);
end

end %main fn


%% Local functions
function [data, studyinfo] = parseSrc(src, options, varargin)
if isstruct(src) && isfield(src,'time') % data struct (single or array)
  data = src; % if length==1,then likely from SimulateModel call
  studyinfo = [];
elseif ischar(src) %string input
  if options.load_all_data_flag % load data
    [data,studyinfo] = ImportData(src, varargin{:});
  else % only load studyinfo
    data = [];
    studyinfo = CheckStudyinfo(src);
  end

  % update study_dir
  if exist(src,'file') && strfind(src, 'studyinfo') %studyinfo.mat
    studyinfo.study_dir = fileparts2(src);
  elseif isdir(src) % study_dir
    studyinfo.study_dir = src;
  end
end

% Old Verbose Way with unnecessary checks
% determine type of src
% if ischar(src)
%   if exist(src,'file') % data file or studyinfo.mat
%     if strfind(src, 'studyinfo') %studyinfo.mat
%       [data,studyinfo] = ImportData(src, varargin{:}); % load data
%       studyinfo.study_dir = fileparts2(src);
%     else % data file
%       [data,studyinfo] = ImportData(src, varargin{:}); % load data
%     end
%   elseif isdir(src) % study_dir
%     [data,studyinfo] = ImportData(src, varargin{:}); % load data
%     studyinfo.study_dir = src;
%   else
%     try
%       [data,studyinfo] = ImportData(src, varargin{:}); % load data
%     catch
%       error('Unknown source for first input/argument.')
%     end
%   end
% elseif isstruct(src) && length(src)>1 % data file cell array
%   data = src;
% elseif isstruct(src) % single data struct or studyinfo struct
%   if isfield(src,'time') % single data file
%     data = src;
%   else % studyinfo struct
%     [data,studyinfo] = ImportData(src, varargin{:}); % load data
%   end
% elseif iscell(src) % cell array of files
%   [data,studyinfo] = ImportData(src, varargin{:}); % load data
% else
%   try
%     [data,studyinfo] = ImportData(src, varargin{:}); % load data
%   catch
%     error('Unknown source for first input/argument.')
%   end
% end
%
%
% % make studyinfo if doesn't exist
% if ~exist('studyinfo','var')
%   studyinfo = [];
% end
end


function [funcIn, nFunc] = parseFuncIn(funcIn)
% convert funcIn input to handle if not a cell array

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
end


function func = parseFunc(func)
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
end


function filename = filenameFromVaried(filename, func, data, plotFnBool, options)
% NOTE: inputs are odd since called from different sources with different
%       states.

if isfield(options, 'save_prefix') && ~isempty(options.save_prefix)
  prefix = options.save_prefix;
else
  if plotFnBool
    if isempty(options.function_options)
      plot_options = options;
    else
      plot_options = options.function_options{fInd};
    end

    % check if 'plot_type' given as part of plot_options
    plot_options = CheckOptions(plot_options,{'plot_type',['plot' num2str(fInd)],[]},false);

    prefix = plot_options.plot_type; % will be waveform by default due to CheckOptions in main fn
  else % ~plotFnBool
    prefix = func2str(func);
  end
end

filename = nameFromVaried(data, prefix, filename);
end


function result = evalFnWithArgs(fInd, data, func, options, varargin)
if isempty(options.function_options)
  parfor dInd = 1:length(data)
    result(dInd) = feval(func,data(dInd),varargin{:});
  end
else
  function_options = options.function_options{fInd};
  parfor dInd = 1:length(data)
    result(dInd) = feval(func,data(dInd),function_options{:});
  end
end
end

function data = loadDataFromSingleSim(simID, options, varargin)
% check if iResult in options.simIDs
if isempty(options.simIDs)
  simIDs = simID;
else
  simIDs = intersect(simID, options.simIDs);

  if isempty(simIDs) % skip if empty
    return
  end
end

varinputs = varargin; % create copy of varargin
if isstruct(varinputs{1})
  varinputs.simIDs = simIDs;
else
  % find simIDs in varinputs
  varinputs{find(~cellfun(@isempty,strfind(varinputs(1:2:end), 'simIDs')))*2} = simIDs;
end

data = ImportData(src, varinputs{:}); % load data
end

function result = add_modifications(result, data)
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
