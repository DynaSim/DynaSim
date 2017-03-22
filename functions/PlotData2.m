function handles=PlotData2(data,varargin)
%% handles=PlotData(data,'option',value)
% Purpose: plot data in various ways depending on what data was provided
% and what options are defined. this function is wrapped by PlotWaveforms,
% PlotPower, ... to provide a single function for organizing and displaying
% data.
% Inputs:
%   data: DynaSim data structure (see CheckData)
%   Accepts the following name/value pairs:
%     'plot_type' {'waveform' (default),'rastergram','rates','power'} - what to plot
%     'population' - name of population to plot (default: 'all'); accepts
%                    regexp strings
%     'variable' - name of variable to plot for each population (default: state variable, X or V);
%                      accepts regexp strings
%     'varied1' - Indices of 1st varied model parameter to plot. If the parameter
%                 is numeric, specify indices (e.g. 1:3 corresponds to 1st-3rd
%                 varied values. If the parameter is a char array, uses
%                 regular expressions. Instead of 'varied1', can also use
%                 the actual parameter name (e.g. 'E_Iapp')
%     'varied2' - As varied 1, for 2nd varied parameter
%     'num_embedded_subplots' - maximum # of waveforms to overlay per plot
%     'max_num_overlaid' - maximum # of waveforms to overlay per plot
%     'do_mean' - {false, true} - Turn on/off averaging across all units
%                 in a population
%     'force_overlay' - {'none', 'populations' (default), 'variables', 'varied1' ... 'variedN'}
%                       If there is only one cell in a population, this forces
%                       PlotData2 to add other information to the overlay.
%     'xlims' - [XMIN XMAX], x-axis limits (default: all data)
%     'ylims' - [YMIN YMAX], y-axis limits (default: all data)
%     'lock_axes' - {false, true}, locks abscissa and ordinate to be
%                                  the same across all subplots
%     'do_zoom' - {false, true} - Turn on zoom function in subplot_grid
%     'yscale' {'linear','log','log10'}, whether to plot linear or log scale
%     'visible' {'on','off'}
%     NOTE: analysis options available depending on plot_type
%       see see CalcFR options for plot_type 'rastergram' or 'rates'
%       see CalcPower options for plot_type 'power'
% Outputs:
%   handles: graphic handles to figures
% 
% See also: CalcFR, CalcPower, PlotWaveforms, CheckData

if ischar(data)
    study_dir = data;
    
    % Import plot files
    data_img = ImportPlots(study_dir);
    
    handles = PlotData2(data_img,varargin{:});
    return;
end

if ~isfield(data,'plot_files')      % Standard DynaSim data structure

    % Check inputs
    data=CheckData(data);
      % note: calling CheckData() at beginning enables analysis/plotting functions to
      % accept data matrix [time x cells] in addition to DynaSim data structure.


    % Convert input data to xPlt
    xp = DynaSim2xPlt(data);
    is_image = 0;
else                            % Structure of links to plots
    
    data_img=data;
    % Load into DynaSim structure
    [data_table,column_titles] = DataField2Table (data_img,'plot_files');

    % The entries in the first column contain the paths to the figure files.
    % There can be multiple figures associated with each simulation, which is
    % why these are cell arrays of strings.
    disp(data_table{1}{1})
    disp(data_table{1}{2})

    % Import the linear data into an xPlt object
    xp = xPlt;
    X = data_table{1}; axislabels = data_table(2:end);
    xp = xp.importLinearData(X, axislabels{:});
    xp = xp.importAxisNames(column_titles(2:end));
    
    % Add dummy population and variable dimensions since the code below
    % expects it
    xd = xp.data;
    xv = xp.exportAxisVals;
    xn = xp.exportAxisNames;
    
    xv(end+1:end+2) = {{'Pop1'},{'X'}};
    xn(end+1:end+2) = {'populations','variables'};
    
    xp = xp.importData(xd,xv);
    xp = xp.importAxisNames(xn);
    clear xd xv xn
    
    % Set up metadata
    % Store metadata info
    meta = struct;
    meta.datainfo(1:2) = nDDictAxis;
    meta.datainfo(1).name = 'time(ms)';
    meta.datainfo(1).values = 1:10;
    meta.datainfo(2).name = 'cells';
        cell_names = [1:5];
        cell_names_str = cellfunu(@(s) ['Cell ' num2str(s)], num2cell(cell_names));
    xp.meta = meta;
    
    
    is_image = 1;

end

% Find out names of varied variables
all_names = xp.exportAxisNames;
varied_names = only_varieds(all_names);  % Returns only the names of the varied variables

% Convert 'varied1'...'variedN' values in varargin to the names of the
% actual varied parameters
myargin = varargin;
for i = 1:length(myargin)
    % Char entries 
    if ischar(myargin{i})
        myargin{i} = variedN_to_axisnames(myargin{i},varied_names);
    end
    
    % Nested char entries within cell array entries
    if iscell(myargin{i})
        for j = 1:length(myargin{i})
            if ischar(myargin{i}{j})
                myargin{i}{j} = variedN_to_axisnames(myargin{i}{j},varied_names);
            end
        end
    end
end
  
% Flag for returning error if the user specifies name/value pairs that are not in the
% CheckOptions list
strict_mode = 1;
  
[options, options_extras0] = CheckOptions(myargin,{...
  'population',[],[],...        
  'variable',[],[],...        
  'num_embedded_subplots',2,{1,2,3,4},...
  'max_num_overlaid',50,[],...
  'do_mean',false,[false true],...
  'force_overlay',[],[],...
  'do_overlay_shift',false,[false true],...
  'overlay_shift_val',[],[],...
  'do_zscore',[false],[false true],...
  'plot_type','waveform',{'waveform','heatmap','rastergram','raster','power','rates'},...
  'xlims',[],[],...
  'ylims',[],[],...
  'zlims',[],[],...
  'lock_axes',1,[false true],...
  'saved_fignum',[1],[],...
  'plot_options',struct,[],...
  'subplot_options',struct,[],...
  'figure_options',struct,[],...
  'do_zoom',false,[false true],...
  'yscale','linear',{'linear','log','log10','log2'},...
  'visible','on',{'on','off'},...
  'save_figures',false,[false true],...
  'save_figname_path',[],[],...
  'supersize_me',false,[false true],...
  'dim_stacking',[],[],...
  },false);
handles=[];

% Options overwrite
if is_image
    options.force_overlay = 'none';
end

% Clause to fix things incase user sets force_overlay to zero!
if isnumeric(options.force_overlay) && ~isempty(options.force_overlay)
    warning('force overlay should be either none or the name the name of an axis {populations, variables, varied1, ... variedN}');
    if options.force_overlay == 0
        options.force_overlay = 'none';
    end
end

% Pull out fields from options struct
plot_type = options.plot_type;
plot_options = options.plot_options;
subplot_options = options.subplot_options;
figure_options = options.figure_options;
num_embedded_subplots = options.num_embedded_subplots;
do_mean = options.do_mean;
force_overlay = options.force_overlay;

% Add default options to structures
% Plot_options
% Used when running xp_matrix_advancedplot3D
    plot_options = struct_addDef(plot_options,'ylims',options.ylims);
    plot_options = struct_addDef(plot_options,'xlims',options.xlims);
    plot_options = struct_addDef(plot_options,'zlims',options.zlims);
% Used when running xp_plotimage
    plot_options = struct_addDef(plot_options,'saved_fignum',options.saved_fignum);

% Subplot_options
subplot_options = struct_addDef(subplot_options,'subplotzoom_enabled',options.do_zoom);

% Figure options
figure_options = struct_addDef(figure_options,'visible',options.visible);
figure_options = struct_addDef(figure_options,'save_figures',options.save_figures);
figure_options = struct_addDef(figure_options,'save_figname_path',options.save_figname_path);
figure_options = struct_addDef(figure_options,'supersize_me',options.supersize_me);






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
    
    % Update cell numbers metadata
    cell_names = [1:max(cellfun(@(x) size(x,2),xp.data(:)))];
    cell_names_str = cellfunu(@(s) ['Cell ' num2str(s)], num2cell(cell_names));
    xp.meta.datainfo(2).values = cell_names_str;
end


% Average across cells if necessary
if do_mean 
    mydata = xp.data;
    mydata = cellfun(@(x) mean(x,2), mydata,'UniformOutput',0);
    xp.data = mydata;
    xp.meta.datainfo(2).values = {'<Cells>'};
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

% User selection for varied parameters
options_extras = options_extras0;
[chosen_varied , options_extras ]= get_chosen_varied(varied_names,options_extras);

% If any options are still leftover, these are extraneous. Report an error
leftover_fields = fieldnames(options_extras);
if ~isempty(leftover_fields) && strict_mode
    error('The following unrecogized name/value pairs were passed in: %s', sprintf('%s ',leftover_fields{:}));
end

% Convert any "all" strings in chosen_varied to colon operators
inds = cellfun(@ischar,chosen_varied);
chosen_varied(inds) = cellfun(@(s) strrep(s,'all',':'),chosen_varied(inds),'UniformOutput',0);

% Select out chosen data
xp2 = xp(chosen_vars,chosen_pop,chosen_varied{:});

% Assign default value to force_overlay
if isempty(force_overlay) 
    if length(xp2.meta.datainfo(2).values) <= 1
        force_overlay = 'populations';
    else
        force_overlay = 'none';
    end
end

if ~strcmpi(force_overlay,'none')
    ax_num = xp2.findaxis(force_overlay);
    if length(ax_num) > 1; error('force_overlay: Ambiguous axis specified. Can only pack at most 1 axis');
    elseif isempty(ax_num)
        error('force_overlay: Cannot find requested axis');
    end
    
    % Save variables associated with this axis
    packed_vars = xp2.axis(ax_num).values;
    packed_name = xp2.axis(ax_num).name;
    
    % Add <average> symbols if necessary to packed_vars
    cellnames = xp2.meta.datainfo(2).values;
    temp = cellfun(@isempty,strfind(cellnames,'<'));    % Check if originals were averages!
    if any(~temp)
        if isnumeric(packed_vars); packed_vars = cellfunu(@(s) [strrep(packed_name,'_',' ') ' ' num2str(s)],num2cell(packed_vars)); end
        packed_vars = cellfunu(@(s) ['<' s '>'], packed_vars);
    end
    
    % Add this new info to datainfo
    % If there is still more than 1 cell present, bump this data down a
    % dimension.
    only_one_cell = length(xp2.meta.datainfo(2).values) == 1;
    if ~only_one_cell
        xp2.meta.datainfo(end+1) = xp2.meta.datainfo(end);
    end
    xp2.meta.datainfo(end).name = packed_name;
    xp2.meta.datainfo(end).values = packed_vars(:)';
    
    
    % Pack the dimension into the first empty dimension
    xp2 = xp2.packDim(force_overlay);
    if only_one_cell
        xp2.data = cellfunu(@squeeze,xp2.data);
    end

end

if options.do_zscore
    mydata = xp2.data;
    for i = 1:numel(mydata)
        mydata{i} = zscore(mydata{i});
    end
    xp2.data = mydata;
end

% Shift the overlay by a certain amount
if options.do_overlay_shift
    mydata = xp2.data;
    for i = 1:numel(mydata)
        mydata{i} = do_shift_lastdim (mydata{i},options.overlay_shift_val);
    end
    xp2.data = mydata;
end

% Squeeze to eliminate superfluous dimensions
xp2 = xp2.squeeze;
Nd = ndims(xp2);

% Set up legend entries
subplot_options.legend1 = setup_legends(xp2);


% Get axis lims
if isempty(plot_options.xlims) && options.lock_axes && ~is_image
    xdat = xp.meta.datainfo(1).values;
    plot_options.xlims = [min(xdat) max(xdat)];
end
if isempty(plot_options.ylims) && options.lock_axes && ~is_image
    switch plot_type
        case 'waveform'
            % Merge all data into one single huge column
            data_all = xp2.data(:);
            data_all = cellfunu(@(x) x(:), data_all);
            data_all = vertcat(data_all{:});
            % Find the max and minima - these are the largest and smallest
            % values we could ever see.
            data_lims = [min(data_all) max(data_all)];
            plot_options.ylims = data_lims;
    end
end

if isempty(plot_options.zlims) && options.lock_axes && ~is_image
    switch plot_type
        case 'heatmap'
            data_all = [xp2.data{:}];
            data_all = data_all(:);
            data_lims = [min(data_all) max(data_all)];
            plot_options.zlims = data_lims;
    end
end

if is_image
    % Is an image
    data_plothandle = @xp_plotimage;
    plot_options.scale = .5;           % Scale of .5 enforces some anti-aliasing
else
    % Is data
    data_plothandle = @xp_matrix_advancedplot3D;
end

% Split available axes into the number of dimensions supported by each
% axis handle
switch num_embedded_subplots
    case 1
        % Ordering of axis handles
        function_handles = {@xp_handles_newfig, @xp_subplot_grid,data_plothandle};   % Specifies the handles of the plotting functions
        dims_per_function_handle = [1,1,1];
        function_args = {{figure_options},{subplot_options},{plot_options}};
        
    case 2
        % Ordering of axis handles
        function_handles = {@xp_handles_newfig, @xp_subplot_grid,data_plothandle};   % Specifies the handles of the plotting functions
        dims_per_function_handle = [1,2,1];
        function_args = {{figure_options},{subplot_options},{plot_options}};
        
    case 3
        % Ordering of axis handles
        function_handles = {@xp_handles_newfig, @xp_subplot_grid, @xp_subplot_grid,data_plothandle};   % Specifies the handles of the plotting functions
        dims_per_function_handle = [1,2,1,1];
        subplot_options2 = subplot_options;
        subplot_options2.legend1 = [];
        subplot_options.display_mode = 1;
        function_args = {{figure_options},{subplot_options2},{subplot_options},{plot_options}};
    case 4
        % Ordering of axis handles
        function_handles = {@xp_handles_newfig, @xp_subplot_grid, @xp_subplot_grid,data_plothandle};   % Specifies the handles of the plotting functions
        dims_per_function_handle = [1,2,2,1];
        subplot_options2 = subplot_options;
        subplot_options2.legend1 = [];
        subplot_options.display_mode = 1;
        function_args = {{figure_options},{subplot_options2},{subplot_options},{plot_options}};
end

% Linearize dimensions of xp2 that are in excess of the total number we can
% plot
maxNplotdims = sum(dims_per_function_handle)-1;
xp2 = reduce_dims(xp2,maxNplotdims);

% Stack up available dimensions based on how much each axis handle can hold
ax_names = [xp2.exportAxisNames, 'data'];
if ~isempty(options.dim_stacking)
    if length(options.dim_stacking) ~= length(ax_names) -1
        error('Incorrect number of dimensions specified. dim_stacking must be some permutation of the following: %s', sprintf('%s ',ax_names{1:end-1}));
    end
    dimensions = get_dimensions({options.dim_stacking{:}, 'data'},dims_per_function_handle);
else
    
    dimensions = get_dimensions(ax_names,dims_per_function_handle);
end

% Remove any excess function handles that aren't needed
available_dims = ~cellfun(@isempty,dimensions);
function_handles = function_handles(available_dims);
dimensions = dimensions(available_dims);
function_args = function_args(available_dims);

% Open new figure if necessary & plot the data
if ~isequal(@xp_handles_newfig, function_handles{1})
    % Cheap hack to force it to create a new figure using our desired
    % parameters for instances when it wouldn't normally call
    % xp_handles_newfig.
    xp3 = xPlt;
    fhandle = @() recursivePlot(xp2,function_handles,dimensions,function_args);
    xp3 = xp3.importData({fhandle});
    xp_handles_newfig(xp3,figure_options);
else
    xp2.recursivePlot(function_handles,dimensions,function_args);
end


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

function [chosen_varied, options_varied ]= get_chosen_varied(varied_names,options_varied)

    
    
    
    % Initialize output
    chosen_varied = repmat({':'},1,length(varied_names));
    
    % Varied name-value pairs entered by user
    varied_NVPs = fieldnames(options_varied);
    
    % See if any of these match actual varied parameters
    for i =  1:length(varied_NVPs)
        ind = find(strcmp(varied_names,varied_NVPs{i}));
        if length(ind) == 1
            chosen_varied{ind} = options_varied.(varied_NVPs{i});
            
            % Optional (remove from options_varied)
            options_varied = rmfield(options_varied,varied_NVPs{i});
        elseif length(ind) > 1
            error('Multiple varied arguments found');
        else
            % Not a varied variable name
        end
        
        
    end
    
end

function str_out = variedN_to_axisnames(str_in,ax_names_varied)

        if strcmp(str_in(1:min(6,end)),'varied')        % User has entered variedX
            % fn is original fieldname (e.g. variedX)
            % fn2 is new field name of varied parameter (e.g. E_Iapp)
            str_out = ax_names_varied{str2num(str_in(7:end))};
        else
            str_out = str_in;
        end
    
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
        dimensions{i} = ax_names(max(1,end-Ndims_curr+1):end);
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

function varied_names = only_varieds(all_names)
    inds = true(1,length(all_names));
    inds(strcmp(all_names,'populations')) = false; 
    inds(strcmp(all_names,'variables')) = false;
    varied_names = all_names(inds);
end

function mydata_out = do_shift_lastdim (mydata,shift)
    sz = size(mydata);
    nd = ndims(mydata);
    
    if isempty(shift)
        % Do adaptive shift
        upscale_factor = 2;
        temp = reshape(mydata,prod(sz(1:nd-1)),sz(nd));
        stdevs = nanstd(temp)*upscale_factor;               % STD ignoring NaNs
        sh = [0, stdevs(1:end-1) + stdevs(2:end)]';
        sh = sh * -1;        % Forces shifts to be downward (same as subplots)
    else
        sh = shift*[0:sz(end)-1]';      % Fixed shift amount
        sh = sh * -1;        % Forces shifts to be downward (same as subplots)
    end
    
    
    sh = permute(sh, [2:nd,1]);
    if length(sz(1:nd-1)) == 1
        sh2 = repmat(sh, sz(1:nd-1),1);     % Special case for scalar input to repmat. When repmat receives a scalar, it repeats BOTH rows and columns instead of just rows
    else
        sh2 = repmat(sh, sz(1:nd-1));
    end
    
    mydata_out = mydata + sh2;
    
end

function leg1 = setup_legends(xp2)
    
    % Pull out all metadata names and values
    for i = 1:length(xp2.meta.datainfo)
        mn{i} = xp2.meta.datainfo(i).name;
        mv{i} = xp2.meta.datainfo(i).values;
    end

    % Convert any numeric entries to cell strings as needed
    for i = 2:length(xp2.meta.datainfo)
        if isnumeric(mv{i})
            % If axis is numeric, as in the case with varied parameters, convert to
            % a cell array of strings
            mv{i} = cellfun(@num2str,num2cell(mv{i}),'UniformOutput',0);

            % Also pre-pend the name of the variable being varied
            for j = 1:length(mv{i})
                mv{i}{j} = [strrep(mn{i},'_',' ') ' ' mv{i}{j}];
            end
        end
    end
    
    if length(mv) == 2
        leg1 = mv{2};
    elseif length(mv) == 3
        %Cartesean product of mv{2} and mv{3}
        k = 0;
        for j = 1:length(mv{3})
            for i = 1:length(mv{2})         % We plot the 2nd dimension 1st; they are grouped together
                k=k+1;
                leg1{k} = [mv{2}{i} ' ' mv{3}{j}];
            end
        end
    else
        error('should not reach');
    end
end
