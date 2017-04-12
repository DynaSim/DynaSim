function result = dsAnalyze(data,func,varargin)
%DSANALYZE - Apply an analysis function to DynaSim data, optionally saving data
%
% Pass a single DynaSim data structure or an array of data structures to a
% user-specified analysis function, add varied info to the results and
% optionally save the output structure.
%
% Usage:
%   result = dsAnalyze(data,func,'option1',value1,...)
%
% Inputs:
%   - data: DynaSim data structure (also accepted: data file name)
%   - func: function handle pointing to analysis function
%   - options: (key/value pairs are passed on to the analysis function)
%     'save_results_flag'   : whether to save result {0 or 1} (default: 0)
%     'result_file'         : where to save result (default: 'result.mat')
%     'format'              : format for saved plots if figures are generated
%                             {'svg','jpg','eps','png'} (default: 'svg')
%     'varied_filename_flag': whether to make filename based on the varied
%                             parameters and type of plot {0 or 1} (default: 0)
%
% Outputs:
%   - result: structure returned by the analysis function
%
% TODO: annotate figures with data set-specific modifications
%
% See also: dsSimulate

% check inputs
options=ds.checkOptions(varargin,{...
  'result_file','result',[],...
  'save_results_flag',0,{0,1},...
  'format','svg',{'svg','jpg','eps','png','fig'},...
  'varied_filename_flag',0,{0,1},...
  'plot_type','waveform',{'waveform','rastergram','raster','power','rates','imagesc','heatmapFR','heatmap_sortedFR','meanFR','meanFRdens'},...
  'save_prefix',[],[],...
  },false);

% change result_file
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
  options.result_file = ds.nameFromVaried(data, prefix, options.result_file);
end

% load data if input is not a DynaSim data structure
if ~(isstruct(data) && isfield(data,'time'))
  data=ds.importData(data,varargin{:}); % load data
end

% confirm function handle
if ~isa(func,'function_handle')
  error('post-processing function must be supplied as a function handle');
end
% confirm single analysis function
if numel(func)>1
  error('Too many function handles were supplied. dsAnalyze only applies a single function to a single data set.');
end
% save data if no output is requested
if nargout<1
  options.save_results_flag=1;
end

% convert data to double precision before analysis
for j=1:length(data)
  for k=1:length(data(j).labels)
    fld=data(j).labels{k};
  data(j).(fld)=double(data(j).(fld));
  end
end

% do analysis
fprintf('\tExecuting post-processing function: %s\n',func2str(func));
tstart=tic;
result=feval(func,data,varargin{:});
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
if all(ishandle(result)) % analysis function returned a graphics handle
  for i=1:length(result)
    % save plot
    if options.save_results_flag
      extension=['.' options.format];%'.svg'; % {.jpg,.svg}
      % temporary default: jpg
      if length(result)==1
        fname=[options.result_file extension];
      else
        fname=[options.result_file '_page' num2str(i) extension];
      end
      set(gcf,'PaperPositionMode','auto');
      fprintf('\t\tSaving plot: %s\n',fname);
      
      switch extension
        case '.svg'
          plot2svg(fname,result(i));
        case '.jpg'
          print(result(i),fname,'-djpeg');
        case '.eps'
          print(result(i),fname,'-depsc');
        case '.png'
          print(result(i),fname,'-dpng');
        case '.fig'
          savefig(result(i),fname);
      end
    end
  end
else % analysis function returned derived data
  if isstruct(result)
    result = add_modifications(result);
    
    for i=1:length(result)
      % add options to result structure
      if length(varargin)>1
        for j=1:2:length(varargin)
          result(i).options.(varargin{j})=varargin{j+1};
        end
      else
        result(i).options=[];
      end
    end
  end
  % save derived data
  if options.save_results_flag
    fname=options.result_file;
    extension = '.mat';
    if ~strcmp(fname(end-3:end), extension) %check for .mat extension
      fname=[fname extension];
    end
    fprintf('\t\tSaving derived data: %s\n', fname);
    save(fname,'result','-v7.3');
  end
end

  function result = add_modifications(result)
    % add modifications to result structure, excluding modifications made
    % within experiments. note: while this nested function is similar to 
    % prepare_varied_metadata in dsSimulate, the data structure contains
    % all modifications (those within and across experiments; listed in 'varied'). 
    % the result structure collapses data sets from an experiment into a single
    % result; thus, each result corresponds to modifications across
    % experiments but not within them; those modifications are stored in
    % the simulator options.
    if ~isempty(data(1).simulator_options.modifications)
      varied={};
      mods=data(1).simulator_options.modifications;
      for ii=1:length(result)
        for jj=1:size(mods,1) 
          % prepare valid field name for thing varied:
          fld=[mods{jj,1} '_' mods{jj,2}];
          
          % convert arrows and periods to underscores
          fld=regexprep(fld,'(->)|(<-)|(-)|(\.)','_');
          
          % remove brackets and parentheses
          fld=regexprep(fld,'[\[\]\(\)\{\}]','');
          result(ii).(fld)=mods{jj,3};
          varied{end+1}=fld;
        end
        result(ii).varied=varied;
        result(ii).modifications=mods;
      end
    elseif isfield(data,'varied') && length(data)==1
      % add 'varied' info from data to result structure
      for ii=1:length(result)
        result(ii).varied=data(1).varied;
        for jj=1:length(data(1).varied)
          result(ii).(data(1).varied{jj})=data(1).(data(1).varied{jj});
        end
      end
    end
  end

end
