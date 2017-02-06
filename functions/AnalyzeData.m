function result=AnalyzeData(data,func,varargin)
%% result=AnalyzeData(data,func,'option1',value1,...)
% purpose: pass a single DynaSim data structure or an array of data
% structures to a user-specified analysis function, add varied info to the
% results and optionally save the output structure.
% 
% inputs:
%   data: DynaSim data structure (also accepted: data file name)
%   func: function handle pointing to analysis function
%   options: 
%     - key/value pairs passed on to the analysis function
%     - 'save_results_flag' (0 or 1) (default: 0): whether to save result
%     - 'result_file' (default: 'result.mat'): where to save result
%     - 'format' (default: 'svg'): format for saved plots if figures are generated
%       {'svg','jpg','eps','png'}
%     - 'varied_filename_flag' (0 or 1) (default: 0): whether to make filename
%       based on the varied parameters and type of plot
% 
% outputs:
%   result: structure returned by the analysis function
% 
% see also: SimulateModel

% todo: annotate figures with data set-specific modifications

% check inputs
options=CheckOptions(varargin,{...
  'result_file','result',[],...
  'save_results_flag',0,{0,1},...
  'format','svg',{'svg','jpg','eps','png'},...
  'varied_filename_flag',0,{0,1},...
  'plot_type','waveform',{'waveform','rastergram','raster','power','rates'},...
  },false);

% change result_file
if options.varied_filename_flag && isfield(data, 'varied')
  options.result_file = nameFromVaried(data, options.plot_type, options.result_file);
end

% load data if input is not a DynaSim data structure
if ~(isstruct(data) && isfield(data,'time'))
  data=ImportData(data,varargin{:}); % load data
end

% confirm function handle
if ~isa(func,'function_handle')
  error('post-processing function must be supplied as a function handle');
end
% confirm single analysis function
if numel(func)>1
  error('Too many function handles were supplied. AnalyzeData only applies a single function to a single data set.');
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
fprintf('Executing post-processing function: %s\n',func2str(func));
tstart=tic;
result=feval(func,data,varargin{:});
fprintf('Elapsed time: %g sec\n',toc(tstart));

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
      fprintf('Saving plot: %s\n',fname);
      switch extension
        case '.svg'
          plot2svg(fname,result(i));
        case '.jpg'
          print(result(i),fname,'-djpeg');
        case '.eps'
          print(result(i),fname,'-depsc');
        case '.png'
          print(result(i),fname,'-dpng');
      end
    end
  end
else % analysis function returned derived data
  if isstruct(result)
    add_modifications;
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
    fprintf('Saving derived data: %s\n',options.result_file);
    save(options.result_file,'result','-v7.3');
  end
end

  function add_modifications
    % add modifications to result structure, excluding modifications made
    % within experiments. note: while this nested function is similar to 
    % prepare_varied_metadata in SimulateModel, the data structure contains
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

