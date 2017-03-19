function handles=PlotData2(data,varargin)
%% handles=PlotData(data,'option',value)
% Purpose: plot data in various ways depending on what data was provided
% and what options are defined. this function is wrapped by PlotWaveforms,
% PlotPower, ... to provide a single function for organizing and displaying
% data.
% Inputs:
%   data: DynaSim data structure (see CheckData)
%   options:
%     'plot_type' {'waveform' (default),'rastergram','rates','power'} - what to plot
%     'variable' - name of field containing data to plot
%                  (default: all pops with state variable of variable in data.labels)
%     'time_limits' - [beg,end] (units of data.time)
%     'max_num_overlaid' - maximum # of waveforms to overlay per plot
%     'max_num_rows' - maximum # of subplot rows per figure
%     'xlim' - [XMIN XMAX], x-axis limits (default: all data)
%     'yscale' {'linear','log','log10'}, whether to plot linear or log scale
%     'visible' {'on','off'}
%     NOTE: analysis options available depending on plot_type
%       see see CalcFR options for plot_type 'rastergram' or 'rates'
%       see CalcPower options for plot_type 'power'
% Outputs:
%   handles: graphic handles to figures
% 
% See also: CalcFR, CalcPower, PlotWaveforms, CheckData

% Check inputs
data=CheckData(data);
  % note: calling CheckData() at beginning enables analysis/plotting functions to
  % accept data matrix [time x cells] in addition to DynaSim data structure.

options=CheckOptions(varargin,{...
  'time_limits',[-inf inf],[],...
  'population',[],[],...        
  'variable',[],[],...        
  'varied',{[]},[],...
  'max_num_overlaid',50,[],...
  'plot_type','waveform',{'waveform','waveform_mean','rastergram','raster','power','rates'},...
  'xlim',[],[],...
  'ylim',[],[],...
  'plot_options',struct,[],...
  'subplot_options',struct,[],...
  'do_zoom',1,[0 1],...
  'yscale','linear',{'linear','log','log10','log2'},...
  'visible','on',{'on','off'},...
  'save_data','on',{'on','off'},...
  },false);
handles=[];

% Pull out fields from options struct
plot_options = options.plot_options;
subplot_options = options.subplot_options;


% Add default options to structures
% Plot_options
plot_options = struct_addDef(plot_options,'ylims',options.ylim);
plot_options = struct_addDef(plot_options,'xlims',options.xlim);

% Subplot_options
subplot_options = struct_addDef(subplot_options,'subplotzoom_enabled',options.do_zoom);


% todo: add option 'plot_mode' {'trace','image'}

time=data.time;

% do any analysis if necessary and set x-data
switch options.plot_type
  case 'waveform'   % plot VARIABLE
    xdata=time;
    xlab='time (ms)'; % x-axis label
  case 'power'      % plot VARIABLE_Power_SUA.Pxx
    if any(cellfun(@isempty,regexp(var_fields,'.*_Power_SUA$')))
      data=CalcPower(data,varargin{:});
    end
    xdata=data(1).([var_fields{1} '_Power_SUA']).frequency;
    xlab='frequency (Hz)'; % x-axis label
    % set default x-limits for power spectrum
    if isempty(options.xlim)
      options.xlim=[0 200]; % Hz
    end
  case {'rastergram','raster'} % raster VARIABLE_spike_times
    if any(cellfun(@isempty,regexp(var_fields,'.*_spike_times$')))
      data=CalcFR(data,varargin{:});
    end
    xdata=time;
    xlab='time (ms)'; % x-axis label
  case 'rates'      % plot VARIABLE_FR
    if any(cellfun(@isempty,regexp(var_fields,'.*_FR$')))
      data=CalcFR(data,varargin{:});
    end
    xdata=data.time_FR;
    xlab='time (ms, bins)'; % x-axis label
end
% if isempty(options.xlim)
%   options.xlim=[min(xdata) max(xdata)];
% end

% Extract the data in a linear table format
[data_table,column_titles,time] = Data2Table (data);

% % Preview the contents of this table
% %     Note: We cannot make this one big cell array since we want to allow
% %     axis labels to be either strings or numerics.
% previewTable(data_table,column_titles);

% Import the linear data into an xPlt object
xp = xPlt;
X = data_table{1};                          % X holds the data that will populate the multidimensional array. Must be numeric or cell array.
axislabels = data_table(2:end);             % Each entry in X has an associated set of axis labels, which will define its location in multidimensional space. **Must be numeric or cell array of chars only**
xp = xp.importLinearData(X,axislabels{:});
xp = xp.importAxisNames(column_titles(2:end));  % There should be 1 axis name for every axis, of type char.


% Apply max overlaid
MTPP = options.max_num_overlaid; % max traces per plot
if any(strcmp(options.plot_type,{'waveform','power'})) && all(cellfun(@isnumeric,xp.data(:)))
    mydata = xp.data;
    mydata2 = cell(size(mydata));
    for i = 1:numel(mydata)
        if ~isempty(mydata{i})
            mydata2{i} = mydata{i}(:,1:min(size(mydata{i},2),MTPP));
        end
    end
    xp.data = mydata2;
    clear mydata mydata2
end


% Axis indices of populations
ax_ind_var = xp.findaxis('variables');
ax_ind_pop = xp.findaxis('population');
ax_ind_varied = get_ax_ind_varied(xp);
ax_ind_varied = find(ax_ind_varied);

% Permute to put varied variables last
xp = permute(xp,[ax_ind_var, ax_ind_pop, ax_ind_varied(:)']);


% User selection for populations
chosen_pop = options.population;
if isempty(chosen_pop)
    chosen_pop = ':';
end

% User selection for variables
chosen_vars = options.variable;
if isempty(chosen_vars)
    chosen_vars = getdefaultvar(xp);
end

% User selection for remaining varied parameters
chosen_varied = get_chosen_varied(xp,options.varied);

% Select out chosen data
xp2 = xp(chosen_vars,chosen_pop,chosen_varied{:});


% % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % Insert code for averaging, etc, here
% % % % % % % % % % % % % % % % % % % % % % % % % 


% Squeeze to eliminate superfluous dimensions
xp2 = xp2.squeeze;
Nd = ndims(xp2);


% If still have too many dimensions, then linearize varied dimensions
% together
maxNplotdims = 3;



 

num_embedded_subplots = 2;

% Split available axes into the number of dimensions supported by each
% axis handle
switch num_embedded_subplots
    case 1
        % Ordering of axis handles
        function_handles = {@xp_handles_newfig, @xp_subplot_grid,@xp_matrix_basicplot};   % Specifies the handles of the plotting functions
        dims_per_function_handle = [1,1,1];
        function_args = {{},{subplot_options},{plot_options}};
        
    case 2
        % Ordering of axis handles
        function_handles = {@xp_handles_newfig, @xp_subplot_grid,@xp_matrix_basicplot};   % Specifies the handles of the plotting functions
        dims_per_function_handle = [1,2,1];
        function_args = {{},{subplot_options},{plot_options}};
        
        
    case 3
        % Ordering of axis handles
        function_handles = {@xp_handles_newfig, @xp_subplot_grid, @xp_subplot_grid,@xp_matrix_basicplot};   % Specifies the handles of the plotting functions
        dims_per_function_handle = [1,2,1,1];
        subplot_options2.display_mode = 1;
        function_args = {{},{},{subplot_options2},{plot_options}};
    case 4
        % Ordering of axis handles
        function_handles = {@xp_handles_newfig, @xp_subplot_grid, @xp_subplot_grid,@xp_matrix_basicplot};   % Specifies the handles of the plotting functions
        dims_per_function_handle = [1,2,2,1];
        subplot_options2.display_mode = 1;
        function_args = {{},{},{subplot_options2},{plot_options}};
end

maxNplotdims = sum(dims_per_function_handle)-1;
xp2 = reduce_dims(xp2,maxNplotdims);

ax_names = [xp2.exportAxisNames, 'data'];
dimensions = get_dimensions(ax_names,dims_per_function_handle);


available_dims = ~cellfun(@isempty,dimensions);
function_handles = function_handles(available_dims);
dimensions = dimensions(available_dims);
function_args = function_args(available_dims);

xp2.recursivePlot(function_handles,dimensions,function_args);


end

function vars_out = getdefaultvar(xp)
    % search through and try to find the variable represnting voltage. If can't find
    % it, just return the first variable listed.
    
    
    % Pull out variables
    vars_orig = xp.axis('variables').values;
    
    % Make everything uppercase to ensure
    % case-insensitive.
    vars = upper(vars_orig);
    possibilities = upper({'V','X','Vm','Xm','Y','Ym'});
    
    ind = [];
    i=0;
    while isempty(ind) && i < length(possibilities)
        i=i+1;
        ind = find(strcmpi(vars,possibilities{i}));
    end
    
    if ~isempty(ind)
        vars_out = vars_orig{ind};
    else
        vars_out = vars_orig{1};
    end
end


function chosen_varied = get_chosen_varied(xp,varied)
    
    chosen_varied = repmat({':'},1,xp.ndims);
    if iscell(varied{1})
        for i = 1:length(varied)
            chosen_varied{xp.findaxis(varied{i}{1})} = varied{i}{2};
        end
    elseif ~isempty(varied{1})
        chosen_varied{xp.findaxis(varied{1})} = varied{2};
    else
    end
    
    chosen_varied = chosen_varied(3:end);       % Drop 1st two entries, since these 
                                                % correspond to population
                                                % and variables, which
                                                % should always be empty.
end


function ax_ind_varied = get_ax_ind_varied(xp)
    % get logical indices of axess corresponding to varied parameters
    ax_ind_pop = xp.findaxis('population');
    ax_ind_var = xp.findaxis('variables');
    ax_ind_varied = true(1,ndims(xp));
    if ~isempty(ax_ind_pop); ax_ind_varied(ax_ind_pop) = 0; end
    if ~isempty(ax_ind_var); ax_ind_varied(ax_ind_var) = 0; end
end

function dimensions = get_dimensions(ax_names,dims_per_function_handle)
    % Split available axes into the number of dimensions supported by each
    % axis handle

    i = length(dims_per_function_handle);
    while i > 0 && ~isempty(ax_names)
        
        Ndims_curr = dims_per_function_handle(i);
        dimensions{i} = ax_names(end-Ndims_curr+1:end);
        ax_names = ax_names(1:end-Ndims_curr);
        i=i-1;
    end
    
end



function xp2 = reduce_dims(xp2,maxNplotdims)
    Nd = ndims(xp2);
    if Nd > maxNplotdims 
        xp2 = xp2.mergeDims( [maxNplotdims:Nd] );
        xp2 = xp2.squeeze;
        Nd = ndims(xp2);

        if Nd ~= maxNplotdims; error('something wrong'); end
    end
end