function result = dsAnalyze(src,varargin)
%DSANALYZE - Apply an analysis function to DynaSim data, optionally saving data
%
% Pass a single DynaSim data structure or an array of data structures to a
% user-specified analysis function, add varied info to the results and
% optionally save the output structure.
%
% Usage:
%  1) dsAnalyze Style: explicit func handle/cell of handles
%   result = dsAnalyze(data,func,'option1',value1,...) % pass data or datafile name
%   result = dsAnalyze(studyinfo,func,'option1',value1,...) % pass studyinfo struct
%   result = dsAnalyze(study_dir,func,'option1',value1,...) % pass study_dir containing studyinfo.mat
%  2) dsSimluate Style: implicit func through options.analysis_functions/options.plot_functions/options.result_functions
%   result = dsAnalyze(data,'option1',value1,...) % pass data or datafile name
%   result = dsAnalyze(studyinfo,'option1',value1,...) % pass studyinfo struct
%   result = dsAnalyze(study_dir,'option1',value1,...) % pass study_dir containing studyinfo.mat
%
% Inputs:
%   - First input/argument:
%     - data: DynaSim data structure or data file path(s)
%     - studyinfo: DynaSim studyinfo structure or path to studyinfo
%     - study_dir: DynaSim study directory containing studyinfo.mat
%   - func: function handle or cell array of function handles pointing to plot
%           or analysis function(s). Should not contain case-insensitive string
%           'plot' unless is a function that returns a figure handle for
%           plotting, in which case it must have 'plot' string.
%   - options: (key/value pairs are passed on to the analysis function)
%     'save_results_flag'   : whether to save result {0 or 1} (default: 0)
%     'matCompatibility_flag': whether to save mat files in compatible mode, vs to prioritize > 2GB VARs {0 or 1} (default: 1)
%     'overwrite_flag': whether to overwrite existing result files {0 or 1} (default: 0)
%     'result_file'         : where to save result (default: 'result.mat')
%     'format'              : format for saved plots if figures are generated
%                             {'svg','jpg','eps','png'} (default: 'svg')
%     'varied_filename_flag': whether to make filename based on the varied
%                             parameters and type of plot {0 or 1}. will overwrite
%                             if multiple plots of same type (use 'save_prefix' to
%                             avoid overwrite in that case) (default: 0)
%     'save_prefix'         : if 'varied_filename_flag'==1, add a string prefix 
%                             to the name (default: '')
%     'close_fig_flag'      : Close/hide figures as they are created. Can be specified
%                             for all figs or for individual fig using plot_options.
%                             Default value is same as 'save_results_flag'.
%     'load_all_data_flag'  : whether to load all the data in studyinfo
%                             at once {0 or 1} (default: 0)
%     'verbose_flag'  : whether to display informative messages/logs (default: 0)
%     'parfor_flag'   : whether to use parfor to run analysis {0 or 1} (default: 0)
%     'simIDs'        : numeric array of simIDs to analyze (default: [])
%     'studyinfo_arg_flag': whether to add 'studyinfo' as a function option {0 or 1} (default: 0)
%     'simID_arg_flag': whether to add 'simID' as a function option {0 or 1} (default: 0)
%     'argout_as_cell': arg output as cell array {0,1} (default: 0)
%    3 ways to specify functions:
%     1)
%     'function_options'    : cell array of option cell arrays {'option1',value1,...}
%                             in which each cell corresponds to the options for
%                             the corresponding function cell. if only passing a
%                             single func, can specificy function options as
%                             key,val list as varargin for dsAnalyze
%     2.1)
%     'analysis_functions'  : cell array of analysis function handles
%     'analysis_options'    : cell array of option cell arrays {'option1',value1,...}
%     'plot_functions'      : cell array of plot function handles
%     'plot_options'        : cell array of option cell arrays {'option1',value1,...}
%     2.2)
%     'result_functions'    : cell array of function handles. Function names should
%                             not contain case-insensitive string 'plot' unless they
%                             are functions that returns a figure handle for plotting,
%                             in which case it must have 'plot' string.
%     'function_options'    : cell array of option cell arrays {'option1',value1,...}
%
%
% Outputs:
%   - result: for single fn, result is struct, cell array, or cell contents returned by the analysis function
%             if postHoc, cell array of the form result{iFunc} containing previous for each fn
%
% TODO:
%   - annotate figures with data set-specific modifications
%   - multiple figure return
%
%
% See also: dsSimulate
%
% Author: Jason Sherfey, PhD <jssherfey@gmail.com>
% Copyright (C) 2016 Jason Sherfey, Boston University, USA
% Edited by Ben Polletta and Erik Roberts


%% General cases:
%   - data struct (likely from SimualteModel call)
%   - data struct array
%   - studyinfo with load_all_data_flag==0
%   - studyinfo with load_all_data_flag==1



% Dev note: calls to this fn (for dev branch on 2/16/18)
% - from dsSimulate:
%     if options.save_data_flag || options.save_results_flag
%       siminfo=studyinfo.simulations(sim_ind);
%
%       for f=1:length(siminfo.result_functions)
%         tmpresult=dsAnalyze(tmpdata,siminfo.result_functions{f},'result_file',siminfo.result_files{f},'save_data_flag',1,'save_results_flag',1,siminfo.result_options{f}{:},'parfor_flag',options.parfor_flag);
% 
%         % since the plots are saved, close all generated figures
%         if all(ishandle(tmpresult))
%           close(tmpresult);
%         end
%       end
%     else
%
%     if ~isempty(options.analysis_functions) && nargoutmain > 2
%       dsAnalyze(tmpdata, options.analysis_functions, 'result_file',[], 'save_data_flag',0, 'save_results_flag',options.save_results_flag, 'function_options',options.analysis_options, 'parfor_flag',options.parfor_flag);
%     end
% 
%     if ~isempty(options.plot_functions)
%       dsAnalyze(tmpdata, options.plot_functions, 'result_file',[], 'save_data_flag',0, 'save_results_flag',options.save_results_flag, 'function_options',options.plot_options, 'parfor_flag',options.parfor_flag);
%     end
% - from cluster job:
%     dsAnalyze(data,siminfo.result_functions{i},'result_file',siminfo.result_files{i},'save_data_flag',1,siminfo.result_options{i}{:});

%% localfn output
if ~nargin
  output = localfunctions; % output var name specific to this fn
  return
end

%% check input type
if (nargin > 1) && (iscell(varargin{1}) || all(isfunction(varargin{1})))
  funcIn = varargin{1};
  varargin(1) = [];
else
  funcIn = [];
end

%% Check inputs
options=dsCheckOptions(varargin,{...
  'in_sim_flag',0,{0,1},...
  'result_file',[],[],...
  'prefix','study',[],...
  'save_results_flag',0,{0,1},...
  'matCompatibility_flag',1,{0,1},...  % whether to save mat files in compatible mode, vs to prioritize > 2GB VARs
  'overwrite_flag',0,{0,1},... % whether to overwrite existing data
  'format','svg',{'svg','jpg','eps','png','fig'},...
  'varied_filename_flag',0,{0,1},...
  'plot_type','waveform',{'waveform','rastergram','raster','power','rates','imagesc','heatmapFR','heatmap_sortedFR','meanFR','meanFRdens','FRpanel'},...
  'save_prefix',[],[],...
  'result_functions',[],[],...
  'function_options',{},[],...
  'simIDs',[],[],...
  'studyinfo_arg_flag',0,{0,1},...
  'simID_arg_flag',0,{0,1},...
  'load_all_data_flag',0,{0,1},...
  'auto_gen_test_data_flag',0,{0,1},...
  'unit_test_flag',0,{0,1},...
  'parfor_flag',0,{0,1},...     % whether to run analysis in parallel (using parfor)
  'verbose_flag',0,{0,1},...
  'analysis_functions',[],[],...
  'analysis_options',[],[],...
  'plot_functions',[],[],...
  'plot_options',[],[],...
  'close_fig_flag',-1,{0,1},... % close figures as they are created. can be specified for all figs or for individual fig using plot_options. -1 means not set by user.
  'argout_as_cell',0,{0,1},... % guarantee output as cell array and leave mising data as empty cells
  'auto_gen_test_data_flag',0,{0,1},...
  },false);

%% auto_gen_test_data_flag argin
if options.auto_gen_test_data_flag
  varargs = varargin;
  varargs{find(strcmp(varargs, 'auto_gen_test_data_flag'))+1} = 0;
  varargs(end+1:end+2) = {'unit_test_flag',1};
  argin = [{src},{funcIn}, varargs]; % specific to this function
end

%% Parse options
if (options.parfor_flag && ~options.load_all_data_flag)
  warning('Since load_all_data_flag==0, setting parfor_flag==0');
  options.parfor_flag = 0;
end

if (options.load_all_data_flag && ~options.parfor_flag)
  dsVprintf(options, 'Since load_all_data_flag==1, recommend setting parfor_flag==1 for speedup. \n');
end

%% Save data if no output is requested.
if nargout < 1 && ~options.in_sim_flag
  options.save_results_flag = 1;
  dsVprintf(options, 'Setting save_results_flag=1 since no nargout.\n')
end

if options.parfor_flag && ~strcmp(reportUI,'matlab')
  disp('For GNU Octave users: Do not expect any speed up by using DynaSim''s ''parfor_flag''. In GNU Octave, parfor loops currently default to regular for loops.');
end

%% Parse src
[data, studyinfo, options] = parseSrc(src, options, varargin{:});
% Data at this point:
%   - 'data' as single struct or struct array, or empty
%   - 'studyinfo' struct or empty

if options.load_all_data_flag
  assert(~isempty(data));
  
  % convert data to double precision before analysis
  if ~options.in_sim_flag
    dsVprintf(options, 'Converting data to double precision before analysis.\n')
  end
  data = convertDoublePrecision(data);
end

if options.studyinfo_arg_flag
  options.studyinfo = studyinfo;
end

% check if studyinfo found
if isempty(studyinfo) % empty inSim
  studyinfoBool = false;
  
  % these must be like this if not studyinfo
  % TODO: add check and warning message
  options.studyinfo_arg_flag = 0;
  options.simID_arg_flag = 0;
  options.simIDs = [];
else
  studyinfoBool = true;
  
  if ~isempty(options.simIDs)
    % filter simulations in studyinfo to only given simIDs
    studyinfo.simulations = studyinfo.simulations(options.simIDs);
  end
  
  options.simIdVec = [studyinfo.simulations.sim_id];
end

% check if study_dir defined
if ~isfield(studyinfo,'study_dir') || isempty(studyinfo.study_dir) || ~isdir(studyinfo.study_dir)
  studyinfo.study_dir = pwd;
end

% Data at this point:
%   - 'data' as single struct or struct array, or empty
%   - 'studyinfo' struct with many flds or just 'study_dir' field

%% Parse fn

if isempty(funcIn)
  % make fn fields into cells
  if ~isempty(options.plot_functions) && ~iscell(options.plot_functions)
    options.plot_functions = {options.plot_functions};
  end

  if ~isempty(options.analysis_functions) && ~iscell(options.analysis_functions)
    options.analysis_functions = {options.analysis_functions};
  end
  
  if ~isempty(options.result_functions) && ~iscell(options.result_functions)
    options.result_functions = {options.result_functions};
  end
end

% make funcIn for dsSimulate Style
if ~isempty(funcIn) % style 1
  plotFnBoolVec = [];
  
  % make sure there is one option cell array per function
  if length(options.function_options) < length(funcIn)
    % extend function_options with blank cells
    options.function_options(length(options.function_options)+1:length(funcIn)) = {{}};
  end
elseif ( ~isempty(options.plot_functions) || ~isempty(options.analysis_functions) ) % style 2.1
  % functions
  funcIn = [options.plot_functions(:); options.analysis_functions(:)];
  plotFnBoolVec = false(size(funcIn));
  plotFnBoolVec(1:length(options.plot_functions)) = true; % specify which fn were plot fn
  
  % options
  if isempty(options.plot_options)
    options.plot_options = {{}};
  end
  if isempty(options.analysis_options)
    options.analysis_options = {{}};
  end
  
  % make sure there is one option cell array per function
  if length(options.plot_options) < length(options.plot_functions)
    % extend plot_options with blank cells
    options.plot_options(length(options.plot_options)+1:length(options.plot_functions)) = {{}};
  end
  if length(options.analysis_options) < length(options.analysis_functions)
    % extend analysis_options with blank cells
    options.analysis_options(length(options.analysis_options)+1:length(options.analysis_functions)) = {{}};
  end
  options.function_options = [options.plot_options(:); options.analysis_options(:)];
elseif ~isempty(options.result_functions) % style 2.2
  % functions
  funcIn = options.result_functions;
  plotFnBoolVec = [];
  
  % options
  if isempty(options.result_options)
    options.result_options = {{}};
  end
  
  % make sure there is one option cell array per function
  if length(options.function_options) < length(options.result_functions)
    % extend function_options with blank cells
    options.function_options(length(options.function_options)+1:length(options.result_functions)) = {{}};
  end
end

% convert func string to handle, or check length of cell array
[funcIn, nFunc] = parseFuncIn(funcIn);

% check if postHoc (ie not from dsSimulate call)
postHocBool = ~options.in_sim_flag;

% allFnResults
if nargout
  allFnResults = cell(nFunc,1);
end

% Handle existing results
[lastPlotIndex, lastAnalysisIndex] = findIndexFromExistingResults();

%% Calc results
plotFnInd = lastPlotIndex;
analysisFnInd = lastAnalysisIndex;

for fInd = 1:nFunc % loop over function inputs
  func = funcIn{fInd};
  
  if ~options.in_sim_flag
    dsVprintf(options, 'Function (%i/%i): %s \n', fInd,nFunc,func2str(func));
  end
  % confirm func is function handle or convert to one if possible
  func = parseFunc(func);

  % check if plot in fn name
  if ~isempty(plotFnBoolVec)
    plotFnBool = plotFnBoolVec(fInd);
  else
    plotFnBool = ~isempty(regexpi(func2str(func), 'plot'));
  end
  
  % plot format
  fopts = dsCheckOptions(options.function_options{fInd}, {'format',[],{'svg','jpg','eps','png','fig'}}, 0);
  if isempty(fopts.format)
    plotFormat = options.format;
  else
    plotFormat = fopts.format;
  end

  % do analysis
  dsVprintf(options, '  Executing post-processing function: %s\n',func2str(func));
  tstart = tic;

  %% Eval func
  if ~isempty(data) % either in sim or posthoc with load_all_data_flag 
    result = evalFnWithArgs(fInd, data, func, options, varargin{:});
    
    if isempty(result)
      fprintf(2, 'Warning: Result empty for function, ''%s''. \n', func2str(func));
    end
  else % posthoc without load_all_data_flag
    result = [];
  end

  % calc nResults
  if ~isempty(result)
    nResults = length(result);
  elseif ~isempty(data)
    nResults = length(data);
  elseif studyinfoBool
    nResults = length(studyinfo.simulations);
  else
    error('Cannot determine number of results');
  end

  dsVprintf(options, '    Elapsed time: %g sec\n',toc(tstart));

  % Dave: Not all plotting functions will return a plot handle. For
  % example, dsPlot2 returns a nested structure of figure, axis, and plot
  % handles. This command updates it.
  if isstruct(result) && isfield(result,'hcurr')
    result = result.hcurr;
  end

  % determine if result is a plot handle or derived data
  if ~iscell(result) && ((~isempty(data) && all(ishandle(result))) || plotFnBool) % analysis function returned a graphics handle or has plot in fn name or given in plot_functions field
    %% Plot Function
    % will save plots, else return main fn since plot already open
    
    plotFnInd = plotFnInd + 1;

    % loop through results. all results may exist or need to be made during loop
    for iResult = 1:nResults
      if ~options.in_sim_flag
        dsVprintf(options, '  Result (%i/%i) \n', iResult,nResults);
      end
      
      extension = ['.' plotFormat]; % '.svg'; % {.jpg,.svg}
      
      if ~postHocBool % in sim
        if ~isempty(options.result_file)
          % ensure extension is extension
          fPath = options.result_file;
          [parentPath, filename, orig_ext] = fileparts(fPath);
          if length(orig_ext) > 4 % extra periods in name
            orig_ext = fPath(end-3:end);
            if ~strcmp(orig_ext, extension)
              fPath = [fPath(1:end-4) extension];
            end
          else
            if ~strcmp(orig_ext, extension)
              fPath = fullfile(parentPath, [filename extension]);
            end
          end
        end
        
        thisResult = result;
      elseif studyinfoBool % posthoc with studyinfo
        simID = studyinfo.simulations(iResult).sim_id;
        
        if isfield(options, 'result_file') && ~isempty(options.result_file)
          fPath = options.result_file;
        else
          if options.load_all_data_flag
            thisData = data(iResult);
            
            thisResult = result(iResult);
          else % load data
            thisData = loadDataFromSingleSim(studyinfo, simID, options, varargin{:});
            
            %skip if no data
            if isempty(thisData)
              dsVprintf(options, '  Skipping simID=%i since no data.\n', simID);
              continue
            end
            
            if options.simID_arg_flag
              options.thisSimID = simID;
            end
            
            % calc result for this data
            thisResult = evalFnWithArgs(fInd, thisData, func, options, varargin{:});
          end % if ~options.load_all_data_flag
          
          % make fPath
          fName = [options.prefix '_sim' num2str(simID) '_plot' num2str(plotFnInd) '_' func2str(func)];
          
          % change result_file if varied_filename_flag
          fopts = dsCheckOptions(options.function_options{fInd}, {'varied_filename_flag',0,{0,1}}, 0);
          if isempty(fopts.varied_filename_flag)
            varied_filename_flag = options.varied_filename_flag;
          else
            varied_filename_flag = fopts.varied_filename_flag;
          end
          
          if varied_filename_flag && isfield(thisData, 'varied')
            fName = filenameFromVaried(fInd, fName, func, thisData, plotFnBool, options, varargin{:});
          end % varied_filename_flag
          
          % make fPath
          fDir = fullfile(studyinfo.study_dir, 'postHocPlots');
          if ~exist(fDir,'dir') && options.save_results_flag
            mkdir(fDir)
          end
          
          fPath = fullfile(fDir,[fName extension]);
        end
      else  % posthoc without studyinfo
        thisResult = result(iResult);
        simID = iResult; % for skipping warning
        
        % make fName
        if isfield(options, 'result_file') && ~isempty(options.result_file)
          fPath = options.result_file;
        else
          fName = [options.prefix '_data' num2str(iResult) '_plot' num2str(iResult) extension];
          fPath = fullfile(fDir,[fName extension]);
        end
        
        % make fDir
        fDir = fullfile(studyinfo.study_dir, 'postHocPlots');
        if ~exist(fDir,'dir') && options.save_results_flag
          mkdir(fDir)
        end
      end % if ~postHocBool

      % Data needed for plotting:
      %   - thisResult
      %   - fPath

      %skip if no result
      if isempty(thisResult)
        if ~postHocBool
          dsVprintf(options, '  Skipping since no result.\n');
        else
          dsVprintf(options, '  Skipping id=%i since no result.\n', simID);
        end

        continue
      end % if isempty(thisResult)

      if options.save_results_flag
        set(thisResult, 'PaperPositionMode','auto');
        dsVprintf(options, '    Saving plot: %s\n',fPath);

        switch extension
          case '.svg'
            plot2svg(fPath, thisResult, [], [], [], [], [], false);
          case '.jpg'
            print(thisResult,fPath,'-djpeg');
          case '.eps'
            print(thisResult,fPath,'-depsc');
          case '.png'
            print(thisResult,fPath,'-dpng');
          case '.fig'
            savefig(thisResult,fPath);
          otherwise
            error('Unknown plot extension. Try again with known extension. See help(dsAnalyze)')
        end
      end %save_results_flag
        
      if (options.save_results_flag && (options.close_fig_flag ~= 0)) || options.close_fig_flag==1
        close(thisResult)
      end

      if ~options.load_all_data_flag
        data = [];
      end

      % store result
      if nargout
        if ~options.argout_as_cell && isstruct(result) && isfield(result,'time')
          % dynasim type structure to store as struct array
          allFnResults{fInd}(iResult) = result;
        else
          if iscell(result) && length(result) == 1
            % if single cell result, store as cell array cell
            allFnResults{fInd}(iResult) = result;
          else
            % if not single cell result, store inside cell array cell
            allFnResults{fInd}(iResult) = {result};
          end
        end
      end % if nargout
    end %nResults
    
  else % analysis function returned derived data
    %% Analysis Function
    if isstruct(result)
      result = add_modifications(result, data, varargin{:});

      for iResult = 1:length(result)
        % add options to result structure
        if length(varargin) > 1
          for j = 1:2:length(varargin)
            result(iResult).options.(varargin{j}) = varargin{j+1};
          end
        else
          result(iResult).options = [];
        end
      end %iResult
    end %isstruct
    
    analysisFnInd = analysisFnInd + 1;

    % switch names in postHoc
    if postHocBool
      allResults = result;
      clear result;
    end

    for iResult = 1:nResults
      if ~options.in_sim_flag
        dsVprintf(options, '  Result (%i/%i) \n', iResult,nResults);
      end
      
      if ~postHocBool % in sim
        fPath = options.result_file;

        % ensure extension is '.mat'
        extension = '.mat';
        [parentPath, filename, orig_ext] = fileparts(fPath);
        if ~strcmp(orig_ext, extension) %check for .mat extension
          fPath = [parentPath filename extension];
        end
      elseif studyinfoBool % posthoc with studyinfo
        simID = studyinfo.simulations(iResult).sim_id;
        
        if options.load_all_data_flag
          thisData = data(iResult);
          
          result = allResults(iResult);
        else % load data
          thisData = loadDataFromSingleSim(studyinfo, simID, options, varargin{:});

          %skip if no data
          if isempty(thisData)
            dsVprintf(options, '  Skipping simID=%i since no data.\n', simID);
            continue
          end

          if options.simID_arg_flag
            options.thisSimID = simID;
          end

          % calc result for this data
          result = evalFnWithArgs(fInd, thisData, func, options, varargin{:});
        end % if options.load_all_data_flag
        
        % make fPath
        if isfield(options, 'result_file') && ~isempty(options.result_file)
          fPath = options.result_file;
        else
          fName = [options.prefix '_sim' num2str(simID) '_analysis' num2str(analysisFnInd) '_' func2str(func) '.mat'];
          
          % change result_file if varied_filename_flag
          fopts = dsCheckOptions(options.function_options{fInd}, {'varied_filename_flag',0,{0,1}}, 0);
          if isempty(fopts.varied_filename_flag)
            varied_filename_flag = options.varied_filename_flag;
          else
            varied_filename_flag = fopts.varied_filename_flag;
          end
          
          if varied_filename_flag && isfield(thisData, 'varied')
            fName = filenameFromVaried(fInd, fName, func, thisData, plotFnBool, options, varargin{:});
          end % varied_filename_flag
          
          % make fDir
          fDir = fullfile(studyinfo.study_dir, 'postHocResults');
          if ~exist(fDir,'dir') && options.save_results_flag
            mkdir(fDir)
          end
          
          fPath = fullfile(fDir,fName);
        end
      else  % posthoc without studyinfo
        result = allResults(iResult);
        simID = iResult; % for skipping warning

        % make fName
        if isfield(options, 'result_file') && ~isempty(options.result_file)
          fPath = options.result_file;
        else
          fName = [options.prefix '_data' num2str(iResult) '_analysis' num2str(analysisFnInd) '_' func2str(func) '.mat'];
          
          % make fDir
          fDir = fullfile(studyinfo.study_dir, 'postHocResults');
          if ~exist(fDir,'dir') && options.save_results_flag
            mkdir(fDir)
          end
          
          fPath = fullfile(fDir,fName);
        end
      end % if ~postHocBool

      %skip if no result
      if isempty(result)
        if ~postHocBool
          dsVprintf(options, '  Skipping since no result.\n');
        else
          dsVprintf(options, '  Skipping id=%i since no result.\n', simID);
        end

        continue
      end % if isempty(result)
      
      if options.save_results_flag  
        dsExportData(result, 'filename',fPath, 'result_flag',1, varargin{:});
      end % save_results_flag
        
      if ~options.load_all_data_flag
        data = [];
      end

      % store result
      if nargout
        if ~options.argout_as_cell && isstruct(result) && isfield(result,'time')
          % dynasim type structure to store as struct array
          allFnResults{fInd}(iResult) = result;
        else
          if iscell(result) && length(result) == 1
            % if single cell result, store as cell array cell
            allFnResults{fInd}(iResult) = result;
          else
            % if not single cell result, store inside cell array cell
            allFnResults{fInd}(iResult) = {result};
          end
        end
      end % if nargout
    end % nResults
  end % ishandle(result)
end % fInd

% get output arg for postHoc
if nargout && postHocBool
  % simplify output if possible
  
  for iFunc = 1:nFunc
    % if only 1 data cell enter cell
    if iscell(allFnResults{iFunc}) && (length(allFnResults{iFunc})==1)
      allFnResults{iFunc} = allFnResults{iFunc}{1};
    end
  end
  
  % if only 1 function, enter cell
  if nFunc == 1
    allFnResults = allFnResults{1};
  end
  
  result = allFnResults; % switch names again
  clear allFnResults
  
  % output of form: result{iFunc}{iResult}
end

%% auto_gen_test_data_flag argout
if options.auto_gen_test_data_flag
  argout = {result}; % specific to this function

  %dsUnitSaveAutoGenTestData(argin, argout); % TODO: check if needs to be saveAutoGenTestDir
end

%% Nested fn -------------------------------------------------------------------
  function [lastPlotIndex, lastAnalysisIndex] = findIndexFromExistingResults()
    % Handle existing results
    % Check to see if results files exist. Deal with overwriting them, or
    % increasing index number. Note: if use same fn name, assume its a new call so
    % increment index. Only overwrite if all functions are the same as original.
    if options.save_results_flag && postHocBool
      if studyinfoBool
        if ~isempty(studyinfo.simulations(1).result_functions)
          oldFns = sort(cellfun(@func2str, studyinfo.simulations(1).result_functions, 'Uni',0));
          
          files = studyinfo.simulations(1).result_files; % studyinfoBool with prior result_functions
        else
          files = lscell(studyinfo.study_dir); % studyinfoBool but no result_functions
        end
      else
        files = lscell(studyinfo.study_dir); % no studyinfoBool
      end
      
      % get filename without extension or parent path
      files = cellfun(@fileNameFileparts, files, 'Uni',0);
      
      % get lastPlotIndex
      plotFiles = regexpi(files, '_plot(\d+)_(.+)', 'tokens');
      plotFiles = [plotFiles{:}];
      if ~isempty(plotFiles)
        plotIndFn = cat(1, plotFiles{:});
        
        plotInds = cellfun(@str2double, plotIndFn(:,1));
        
        lastPlotIndex = max(plotInds);
      else
        lastPlotIndex = 0;
      end
      
      % get plots lastPlotIndex
      plotsDir = fullfile(studyinfo.study_dir, 'plots');
      if exist(plotsDir ,'dir')
        plotFiles = lscell(plotsDir);
        plotFiles = regexpi(plotFiles, '_plot(\d+)_(.+)', 'tokens');
        plotFiles = [plotFiles{:}];
        if ~isempty(plotFiles)
          plotIndFn = cat(1, plotFiles{:});
          
          plotInds = cellfun(@str2double, plotIndFn(:,1));
          
          lastPlotIndex = max(lastPlotIndex, max(plotInds));
        else
          lastPlotIndex = lastPlotIndex;
        end
      end
      
      % get postHocPlots lastPlotIndex
      postHocPlotsDir = fullfile(studyinfo.study_dir, 'postHocPlots');
      if exist(postHocPlotsDir ,'dir')
        plotFiles = lscell(postHocPlotsDir);
        plotFiles = regexpi(plotFiles, '_plot(\d+)_(.+)', 'tokens');
        plotFiles = [plotFiles{:}];
        if ~isempty(plotFiles)
          plotIndFn = cat(1, plotFiles{:});
          
          plotInds = cellfun(@str2double, plotIndFn(:,1));
          
          lastPlotIndex = max(lastPlotIndex, max(plotInds));
        else
          lastPlotIndex = lastPlotIndex;
        end
      end
      
      % get lastAnalysisIndex
      analysisFiles = regexpi(files, '_analysis(\d+)_(.+)', 'tokens');
      analysisFiles = [analysisFiles{:}];
      if ~isempty(analysisFiles)
        analysisIndFn = cat(1, analysisFiles{:});
        
        analysisInds = cellfun(@str2double, analysisIndFn(:,1));
        
        lastAnalysisIndex = max(analysisInds);
      else
        lastAnalysisIndex = 0;
      end
      
      % get results lastAnalysisIndex
      resultsDir = fullfile(studyinfo.study_dir, 'results');
      if exist(resultsDir ,'dir')
        analysisFiles = lscell(resultsDir);
        analysisFiles = regexpi(analysisFiles, '_analysis(\d+)_(.+)', 'tokens');
        analysisFiles = [analysisFiles{:}];
        if ~isempty(analysisFiles)
          analysisIndFn = cat(1, analysisFiles{:});
          
          analysisInds = cellfun(@str2double, analysisIndFn(:,1));
          
          lastAnalysisIndex = max(lastAnalysisIndex, max(analysisInds));
        else
          lastAnalysisIndex = lastAnalysisIndex;
        end
      end
      
      % get postHocResults lastAnalysisIndex
      postHocResultsDir = fullfile(studyinfo.study_dir, 'postHocResults');
      if exist(postHocResultsDir ,'dir')
        analysisFiles = lscell(postHocResultsDir);
        analysisFiles = regexpi(analysisFiles, '_analysis(\d+)_(.+)', 'tokens');
        analysisFiles = [analysisFiles{:}];
        if ~isempty(analysisFiles)
          analysisIndFn = cat(1, analysisFiles{:});
          
          analysisInds = cellfun(@str2double, analysisIndFn(:,1));
          
          lastAnalysisIndex = max(lastAnalysisIndex, max(analysisInds));
        else
          lastAnalysisIndex = lastAnalysisIndex;
        end
      end
      
      if ~studyinfoBool
        oldFns = sort([analysisIndFn(:,2); plotIndFn(:,2)]);
      end
      
      if options.overwrite_flag && (~isempty(plotFiles) || ~isempty(analysisFiles))
        % check if all functions the same as old ones
        newFns = sort(cellfun(@func2str, funcIn, 'Uni',0));
        
        if isempty(setdiff(oldFns, oldFns))
          dsVprintf(options, 'Overwriting old results and plot function indicies in new folders starting at index 0.');
          
          lastPlotIndex = 0;
          lastAnalysisIndex = 0;
        else
          dsVprintf(options, 'Not overwriting old results and plot function indicies in new folders since functions are not the same, so incrementing index.');
        end
      end
    else % ~postHocBool
      % setting to avoid errors
      lastPlotIndex = 0;
      lastAnalysisIndex = 0;
    end % if options.save_results_flag && postHocBool
  end % findIndexFromExistingResults
% End Nested Fn ----------------------------------------------------------------

end %main fn




%% Local functions

function [data, studyinfo, options] = parseSrc(src, options, varargin)

% if src is:
% - data: data = src, studyinfo = []
% - studydir and studyinfo: just studyinfo, and if load_all_data_flag, loads data

%% auto_gen_test_data_flag argin
warning('off','catstruct:DuplicatesFound');
options = catstruct(options, dsCheckOptions(varargin,{'auto_gen_test_data_flag',0,{0,1}},false));
warning('on','catstruct:DuplicatesFound');
if options.auto_gen_test_data_flag
  varargs = varargin;
  varargs{find(strcmp(varargs, 'auto_gen_test_data_flag'))+1} = 0;
  varargs(end+1:end+2) = {'unit_test_flag',1};
  argin = [{src},{options}, varargs]; % specific to this function
end

if isempty(src)
  src = pwd;
end


if isstruct(src) && isfield(src,'time') % data struct (single or array)
  data = src; % if length==1,then likely from SimulateModel call
  studyinfo = [];
  
  options.load_all_data_flag = 1; % data has been loaded
elseif ischar(src) %string input
  if options.load_all_data_flag % load data
    [data,~,dataExistBoolVec] = dsImport(src, varargin{:});
    % if any data missing, will return struct with fewer entries, but gives
    % dataExistBoolVec showing which data did exist
    
    studyinfo = dsCheckStudyinfo(src);
    
    if ~all(dataExistBoolVec)
      error('Some data missing and handling this has not been implemented yet. Temporary solution is to pass in simIDs of existing data.')
    end
  else % only load studyinfo
    data = [];
    studyinfo = dsCheckStudyinfo(src);
  end

  % update study_dir
  if exist(src, 'file') && contains(src, 'studyinfo') %studyinfo.mat
    studyinfo.study_dir = fileparts2(src);
  elseif isdir(src) % study_dir
    studyinfo.study_dir = src;
  end
end


%% auto_gen_test_data_flag argout
if options.auto_gen_test_data_flag
  argout = {data, studyinfo}; % specific to this function

  %dsUnitSaveAutoGenTestDataLocalFn(argin, argout); % localfn
end

end % parseSrc


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
  
  % convert to fn handle if strings
  if all(cellfun(@ischar, funcIn))
    funcIn = cellfun(@str2func, funcIn, 'Uni',0);
  end
end
end % parseFuncIn


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
end % parseFunc


function filename = filenameFromVaried(fInd, filename, func, data, plotFnBool, options, varargin)
% NOTE: inputs are odd since called from different sources with different
%       states.

%% auto_gen_test_data_flag argin
warning('off','catstruct:DuplicatesFound');
options = catstruct(options, dsCheckOptions(varargin,{'auto_gen_test_data_flag',0,{0,1}},false));
warning('on','catstruct:DuplicatesFound');
if options.auto_gen_test_data_flag
  varargs = varargin;
  varargs{find(strcmp(varargs, 'auto_gen_test_data_flag'))+1} = 0;
  varargs(end+1:end+2) = {'unit_test_flag',1};
  argin = [{filename}, {func}, {data}, {plotFnBool}, {options}, varargs]; % specific to this function
end

fopts = dsCheckOptions(options.function_options{fInd}, {'save_prefix',[],[]}, 0);
if isempty(fopts.save_prefix)
  save_prefix = options.save_prefix;
else
  save_prefix = fopts.save_prefix;
end

if ~isempty(save_prefix)
  prefix = save_prefix;
else
  if plotFnBool
    plot_options = options.function_options{fInd};
    
    if isempty(plot_options)
      plot_options = options;
    end

    % check if 'plot_type' given as part of plot_options
    plot_options = dsCheckOptions(plot_options,{'plot_type',[options.prefix '_' func2str(func)],[]},false);

    prefix = plot_options.plot_type;
  else % ~plotFnBool
    prefix = [options.prefix '_' func2str(func)];
  end
end

filename = dsNameFromVaried(data, prefix, filename);

%% auto_gen_test_data_flag argout
if options.auto_gen_test_data_flag
  argout = {filename}; % specific to this function

  %dsUnitSaveAutoGenTestDataLocalFn(argin, argout); % localfn
end

end % filenameFromVaried



function result = evalFnWithArgs(fInd, data, func, options, varargin)
% if not load_all_data_flag, will be only 1 dataset

try
  make_invis_bool = options.save_results_flag && (options.close_fig_flag ~= 0);
  
  if options.studyinfo_arg_flag
    options.function_options{fInd}(end+1:end+2) = {'studyinfo',options.studyinfo};
  end
  
  if options.simID_arg_flag
    if options.load_all_data_flag
      options.function_options{fInd}(end+1) = {'simID'};
    else
      options.function_options{fInd}(end+1:end+2) = {'simID',options.thisSimID};
    end
  end
  
  if isempty(options.function_options) || isempty(options.function_options{fInd})
    % Only do parfor mode if parfor_flag is set and parpool is already running. Otherwise, this will add unnecessary overhead.
    if options.parfor_flag % && ~isempty(p)
      parfor iData = 1:length(data)
        result(iData) = feval(func,data(iData),varargin{:});
        
        if ishandle(result(iData)) && make_invis_bool
          set(result(iData), 'Visible', 'off'); % cannot close yet until save, but can make invisible
        end
      end
    else
      for iData = 1:length(data)
        result(iData) = feval(func,data(iData),varargin{:});
        
        try
            % For dsPlot2
            if ishandle(result(iData).hcurr) && make_invis_bool
              set(result(iData).hcurr, 'Visible', 'off'); % cannot close yet until save, but can make invisible
            end
        catch
            % For dsPlot
            if ishandle(result(iData)) && make_invis_bool
              set(result(iData), 'Visible', 'off'); % cannot close yet until save, but can make invisible
            end
        end
      end
    end % options.parfor_flag && ~isempty(p)
  else % ~isempty(options.function_options)
    function_options = options.function_options{fInd};
    
    if options.load_all_data_flag && options.simID_arg_flag
      addSimIdBool = 1;
      
      simIdVec = options.simIdVec;
    else
      addSimIdBool = 0;
      
      simIdVec = []; % define so parfor doesnt throw warning
    end

    if options.parfor_flag % && ~isempty(p)
      parfor iData = 1:length(data)
        if ~addSimIdBool
          result(iData) = feval(func, data(iData), function_options{:}); %#ok<PFBNS>
        else
          result(iData) = feval(func, data(iData), function_options{:}, simIdVec(iData));
        end
        
        if ishandle(result(iData)) && make_invis_bool
          set(result(iData), 'Visible', 'off'); % cannot close yet until save, but can make invisible
        end
      end
    else
      for iData = 1:length(data)
        if ~addSimIdBool
          result(iData) = feval(func, data(iData), function_options{:});
        else
          result(iData) = feval(func, data(iData), function_options{:}, simIdVec(iData));
        end
        
        if ishandle(result(iData)) && make_invis_bool
          set(result(iData), 'Visible', 'off'); % cannot close yet until save, but can make invisible
        end
      end
    end % options.parfor_flag && ~isempty(p)
  end % isempty(options.function_options)
catch err
  warning(err.message);
  result = [];
end

end % evalFnWithArgs


function data = loadDataFromSingleSim(src, simID, options, varargin)

% Dev Note: this block is probably not needed anymore, but should not break anything
% check if iResult in options.simIDs
if isempty(options.simIDs)
  simIDs = simID;
else
  simIDs = intersect(simID, options.simIDs);

  if isempty(simIDs) % skip if simID not in options.simIDs
    return
  end
end

% check if src is data or other
if isstruct(src) && isfield('time','src')
  data = src;
  data = data(simID); % this may not work if a subset of data is given
else
  % overwrite simIDs in varargin
  varargin(end+1:end+2) = {'simIDs', simIDs};

  data = dsImport(src, varargin{:}); % load data
end

data = convertDoublePrecision(data);
end % loadDataFromSingleSim


function result = add_modifications(result, data, varargin)
% add modifications to result structure, excluding modifications made
% within experiments. note: while this nested function is similar to
% dsModifications2Vary called by dsSimulate, the data structure contains
% all modifications (those within and across experiments; listed in 'varied').
% the result structure collapses data sets from an experiment into a single
% result; thus, each result corresponds to modifications across
% experiments but not within them; those modifications are stored in
% the simulator options.

%% auto_gen_test_data_flag argin
options = dsCheckOptions(varargin,{'auto_gen_test_data_flag',0,{0,1}},false);
if options.auto_gen_test_data_flag
  varargs = varargin;
  varargs{find(strcmp(varargs, 'auto_gen_test_data_flag'))+1} = 0;
  varargs(end+1:end+2) = {'unit_test_flag',1};
  argin = [{result}, {data}, varargs]; % specific to this function
end

% #todo: The function dsModifications2Vary implements this functionality.
% Consider using it here.
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

%% auto_gen_test_data_flag argout
if options.auto_gen_test_data_flag
  argout = {result}; % specific to this function

  %dsUnitSaveAutoGenTestDataLocalFn(argin, argout); % localfn
end
end % add_modifications


function filename = fileNameFileparts(filepath)
  [~, filename] = fileparts(filepath);
end % fileNameFileparts


function data = convertDoublePrecision(data)
  if ~isempty(data) && isstruct(data)
    for j = 1:length(data)
      for k = 1:length(data(j).labels)
        fld = data(j).labels{k};
        data(j).(fld) = double(data(j).(fld));
      end
    end
  end
end % convertDoublePrecision
