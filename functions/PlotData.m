function handles=PlotData(data,varargin)
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
% Plots:
% if Nsims>1: one sim per row
% elseif Npops>1: one pop per row
% else: one cell per row
% 
% Examples for specifying 'variable' option: 
% []      : plot all data.labels with same variable name as first element of data.labels (eg, 'pop1_v' and 'pop2_v')
% '*'     : plot all data.labels
% '*_v'   : plot all data.labels ending in _v (i.e., all state variables 'v' for all populations)
% 'pop1_*': plot all data.labels starting with pop1_ (i.e., all variables for population 'pop1')
% 'pop1_v': plot only variable 'pop1_v'
% 'v'     : look for all data.labels ending in _v then starting with v_ (eg, all populations with variable 'v')
% 'pop1'  : look for all data.labels ending in _pop1 then starting with pop1_ (eg, all variables for population 'pop1')
% '*_iNa_*': plot all data.labels for the 'iNa' mechanism (for all populations)
%   
% Examples:
% One cell:
% data=SimulateModel('dv/dt=@current+10; {iNa,iK}','tspan',[0 500]);
% PlotData(data); % plot first state variable ('v')
% PlotData(data,'variable','*'); % plot all state variables
% PlotData(data,'variable','*','time_limits',[30 60]); % plot all variables and time 30-60ms
% % plot power spectrum
% PlotData(data,'variable','v','plot_type','power');
% PlotData(data,'variable','*','plot_type','power');
% 
% One population: with noisy input
% data=SimulateModel('dv[5]/dt=@current+10*(1+randn(1,Npop)); {iNa,iK}','tspan',[0 250]);
% PlotData(data);
% PlotData(data,'variable','*'); % plot all state variables (all cells)
% PlotData(data,'variable','m'); % plot state variable 'm' (all cells)
% PlotData(data,'variable','*','time_limits',[30 60]); % plot all variables and time 30-60ms
% % plot power spectrum
% PlotData(data,'variable','v','plot_type','power');
% PlotData(data,'variable','*','plot_type','power');
% % plot rastergram
% PlotData(data,'variable','v','plot_type','rastergram');
% PlotData(data,'variable','*','plot_type','rastergram');
% 
% One population varying one parameter (input amplitude):
% eqns='dv[5]/dt=@current+amp*(1+randn(1,Npop)); {iNa,iK}';
% vary={'','amp',[0 10 20]};
% data=SimulateModel(eqns,'vary',vary,'tspan',[0 200]);
% PlotData(data);
% PlotData(data,'variable','m');
% PlotData(data,'variable','*');
% % plot power spectrum
% PlotData(data,'variable','v','plot_type','power');
% % plot rastergram
% PlotData(data,'variable','v','plot_type','rastergram');
% 
% One population varying two parameters (input amplitude and membrane capacitance):
% eqns='dv[5]/dt=@current/Cm+amp*(1+randn(1,Npop)); {iNa,iK}';
% vary={'','Cm',[1 2]; '','amp',[0 10 20]};
% data=SimulateModel(eqns,'vary',vary,'tspan',[0 200]);
% PlotData(data);
% PlotData(data,'variable','*');
% % plot power spectrum
% PlotData(data,'variable','v','plot_type','power');
% % plot rastergram
% PlotData(data,'variable','v','plot_type','rastergram');
% 
% Two populations: noisy input to E and excitatory connection from E to I
% spec=[];
% spec.populations(1).name='E1';
% spec.populations(1).equations='dv[5]/dt=@current+amp*(1+randn(1,Npop)); amp=10; {iNa,iK}';
% spec.populations(2).name='E2';
% spec.populations(2).equations='dv[2]/dt=@current; {iNa,iK}';
% spec.connections(1).direction='E1->E2';
% spec.connections(1).mechanism_list='iAMPA';
% data=SimulateModel(spec,'tspan',[0 200]);
% PlotData(data); % plot first state variable
% PlotData(data,'variable','*'); 
% % plot monitored synaptic current with post-synaptic voltages:
% PlotData(data,'variable',{'E2_v','ISYN'}); 
% % plot monitored synaptic current with pre- and post-synaptic voltages:
% PlotData(data,'variable',{'v','ISYN'}); 
% % plot power spectrum
% PlotData(data,'variable','v','plot_type','power');
% PlotData(data,'variable',{'E2_v','ISYN'},'plot_type','power');
% % plot rastergram
% PlotData(data,'variable','v','plot_type','rastergram');
% 
% Two populations varying one parameter (input amplitude):
% vary={'E1','amp',[0 10 20]};
% data=SimulateModel(spec,'vary',vary,'tspan',[0 200]);
% PlotData(data);
% PlotData(data,'variable','*');
% PlotData(data,'variable','*_iNa_*');
% % plot power spectrum
% PlotData(data,'variable','v','plot_type','power');
% % plot rastergram
% PlotData(data,'variable','v','plot_type','rastergram');
% 
% Two populations varying two parameters (input amplitude and synaptic conductance):
% vary={'E1','amp',[0 10 20]; 'E1->E2','gSYN',[0 .05 .1]};
% data=SimulateModel(spec,'vary',vary,'tspan',[0 200]);
% % plot voltage waveforms
% PlotData(data,'variable','v','plot_type','power');
% % plot voltage power spectrum
% PlotData(data,'variable','v','plot_type','waveform');
% % plot voltage-derived rastergram
% PlotData(data,'variable','v','plot_type','rastergram');
% % more plots
% PlotData(data,'variable','ISYN');
% PlotData(data,'variable','E1_v');
% PlotData(data,'variable','*');
% 
% See also: CalcFR, CalcPower, PlotWaveforms, CheckData

% Check inputs
data=CheckData(data);
  % note: calling CheckData() at beginning enables analysis/plotting functions to
  % accept data matrix [time x cells] in addition to DynaSim data structure.

fields=fieldnames(data);

vary_labels = data(1).varied; % data(1).simulator_options.vary;

no_vary_labels = length(vary_labels);

vary_params = nan(length(data), no_vary_labels);

vary_vectors = cell(no_vary_labels, 1);

vary_lengths = nan(no_vary_labels, 1);

for v = 1:no_vary_labels
   
    vary_params(:, v) = [data.(vary_labels{v})];
    
    vary_vectors{v} = unique(vary_params(:, v));
    
    vary_lengths(v) = length(vary_vectors{v});
    
end

[effective_vary_lengths, ~] = CheckCovary(vary_lengths, vary_params);

dimensions_varied = sum(effective_vary_lengths > 1);

if dimensions_varied > 2
    
    no_figures = prod(effective_vary_lengths(3:end));
    
    data = reshape(data, prod(effective_vary_lengths(3:end)), prod(effective_vary_lengths([1 2])));
    
    figure_params = nan(no_vary_labels - 2, 1);
    
    vary_lengths_cp = cumprod(vary_lengths);
    
    for f = 1:no_figures
        
        figure_params(1) = vary_vectors{3}(mod(f - 1, vary_lengths(3)) + 1);
        
        for v = 4:no_vary_labels
            
            figure_params(v - 2) = vary_vectors{v}(ceil(f/vary_lengths_cp(v - 3)));
            
        end
        
        vary_title = '';
        
        for v = 1:(no_vary_labels - 2)
            
            vary_title = [vary_title, sprintf('%s = %f ', vary_labels{v + 2}, figure_params(v))];
            
        end
    
        handles = PlotData(data(f, :), varargin{:});
        
        for h = 1:length(handles)
        
            % figure(handles(h))
    
            mtit(handles(h), vary_title, 'FontSize', 14, 'yoff', .2)
       
        end
        
    end
    
    return
    
end

options=CheckOptions(varargin,{...
  'time_limits',[-inf inf],[],...
  'variable',[],[],...        
  'max_num_overlaid',50,[],...
  'max_num_rows',20,[],...
  'plot_mode','trace',{'trace','image'},...
  'plot_type','waveform',{'waveform','rastergram','raster','power','rates'},...
  'xlim',[],[],...
  'ylim',[],[],...
  'yscale','linear',{'linear','log','log10','log2'},...
  'visible','on',{'on','off'},...
  },false);
data=CheckData(data);
handles=[];

% todo: add option 'plot_mode' {'trace','image'}

% variables to plot
var_fields=SelectVariables(data(1).labels,options.variable);
tmp=regexp(var_fields,'_(.+)$','tokens','once');
variables=unique([tmp{:}]);

% populations to plot
pop_names={data(1).model.specification.populations.name}; % list of populations
% restrict to populations with variables to plot
pop_indices=[];     % indices of populations to plot
pop_var_indices={}; % indices of var_fields to plot per population
  % pop_var_indices{pop}(variable): index into var_fields
for i=1:length(pop_names)
  % do any variables start with this population name?
  var_inds=find(~cellfun(@isempty,regexp(var_fields,['^' pop_names{i}])));
  if any(var_inds)
    inds=cellfun(@(x)find(~cellfun(@isempty,regexp(var_fields(var_inds),['_' x '$'],'once'))),variables,'uni',0);
    varsel=~cellfun(@isempty,inds);
    fldind=unique([inds{:}]);    
    pop_indices(end+1)=i;
    pop_var_indices{end+1}=nan(1,length(variables));
    pop_var_indices{end}(varsel)=var_inds(fldind);
  end
end
pop_names=pop_names(pop_indices);

% data set info
time=data(1).time; % time vector
pop_sizes=[data(1).model.specification.populations(pop_indices).size];
num_pops=length(pop_names); % number of populations to plot
num_sims=length(data); % number of simulations
num_vars=length(variables);
num_labels=length(var_fields); % number of labels to plot
num_times=length(time);
% set x-axis limits
if isempty(options.time_limits)
  options.time_limits=[min(time) max(time)];
end

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
if isempty(options.xlim)
  options.xlim=[min(xdata) max(xdata)];
end

MRPF = options.max_num_rows; % max rows per fig
MTPP = options.max_num_overlaid; % max traces per plot

% how many plots:
if num_sims==1 && num_pops==1 && num_vars==1
  num_fig_sets=1; num_figs=ceil(pop_sizes/MRPF); num_rows=min(pop_sizes,MRPF);
elseif num_sims==1 && num_pops==1 && num_vars>1
  num_fig_sets=1; num_figs=ceil(num_vars/MRPF); num_rows=min(num_vars,MRPF);
elseif num_sims==1 && num_pops>1 && num_vars==1
  num_fig_sets=1; num_figs=ceil(num_pops/MRPF); num_rows=min(num_pops,MRPF);
elseif num_sims==1 && num_pops>1 && num_vars>1
  num_fig_sets=num_vars; num_figs=ceil(num_pops/MRPF); num_rows=min(num_pops,MRPF);
elseif num_sims>1 && num_pops==1 && num_vars==1
  num_fig_sets=1; num_figs=ceil(num_sims/MRPF); num_rows=min(num_sims,MRPF);
elseif num_sims>1 && num_pops==1 && num_vars>1
  num_fig_sets=num_vars; num_figs=ceil(num_sims/MRPF); num_rows=min(num_sims,MRPF);
elseif num_sims>1 && num_pops>1 && num_vars==1
  num_fig_sets=1; num_figs=ceil(num_sims/MRPF); num_rows=min(num_sims,MRPF);
elseif num_sims>1 && num_pops>1 && num_vars>1
  num_fig_sets=num_vars; num_figs=ceil(num_sims/MRPF); num_rows=min(num_sims,MRPF);
else
  error('unrecognized dimensions');
end

% make subplot adjustments for varied parameters
if num_sims>1 && isfield(data,'varied')
  % collect info on parameters varied
  varied=data(1).varied;
  num_varied=length(varied); % number of model components varied across simulations
  num_sims=length(data); % number of data sets (one per simulation)
  % collect info on parameters varied
  param_mat=zeros(num_sims,num_varied); % values for each simulation
  param_cell=cell(1,num_varied); % unique values for each parameter
  % loop over varied components and collect values
  for j=1:num_varied
    if isnumeric(data(1).(varied{j}))
      param_mat(:,j)=[data.(varied{j})]; % values for each simulation
      param_cell{j}=unique([data.(varied{j})]); % unique values for each parameter
    else
      % todo: handle sims varying non-numeric model components 
      % (eg, mechanisms) (also in PlotFR and SelectData)
    end
  end
  param_size=cellfun(@length,param_cell); % number of unique values for each parameter
  % varied parameter with most elements goes along the rows (everything else goes along columns)
  row_param_index=find(param_size==max(param_size),1,'first');
  row_param_name=varied{row_param_index};
  row_param_values=param_cell{row_param_index};
  num_rows=length(row_param_values);
  %num_cols=num_sims/num_rows;
  num_figs=ceil(num_rows/MRPF);
  % collect sims for each value of the row parameter
  indices={};
  for row=1:num_rows
    indices{row}=find(param_mat(:,row_param_index)==row_param_values(row));
  end
  num_per_row=cellfun(@length,indices);
  num_cols=max(num_per_row);
  sim_indices=nan(num_cols,num_rows);
  % arrange sim indices for each row in a matrix
  for row=1:num_rows
    sim_indices(1:num_per_row(row),row)=indices{row};
  end
%   sim_indices=[];
%   for row=1:num_rows
%     sim_indices=[sim_indices find(param_mat(:,row_param_index)==row_param_values(row))];
%   end
else
  sim_indices=ones(1,num_rows); % index into data array
  num_cols=1;
end

max_legend_entries=10;
for figset=1:num_fig_sets
  for fig=1:num_figs
    ylims=[nan nan];
    % create figure
    handles(end+1)=figure('units','normalized','outerposition',[0 0 1 1],'visible',options.visible);
    % position axes
    haxes=tight_subplot(num_rows,num_cols,[.01 .03],[.05 .01],[.03 .01]);
    axis_counter=0;
    AuxData=[];
    vlines=[];
    allspikes={};
    text_string=''; % string to add to subplot (set below)
    legend_strings=''; % legend for subplot (set below)
    shared_ylims_flag=1;
    % draw plots
    for row=1:num_rows
      for col=1:num_cols
        dat=[];
        sim_index=sim_indices(col,row); % index into data array for this subplot
        axis_counter=axis_counter+1; % number subplot axis we're on
        if isnan(sim_index)
          continue;
        end
        % #################################################################
        % what to plot
        % -----------------------------------------------------------------
        if num_sims==1 && num_pops==1 && num_vars==1
        % -----------------------------------------------------------------
          % one cell per row: dat = data(s=1).(var)(:,c=r) where var=vars{v=1}
          var=var_fields{1};
          switch options.plot_type
            case 'waveform'
              dat=data(sim_index).(var)(:,row);
            case 'power'
              AuxData=data(sim_index).([var '_Power_MUA']).Pxx;
              vlines=data(sim_index).([var '_Power_MUA']).PeakFreq;
              AuxDataName={'MUA Power'};
              var=[var '_Power_SUA'];
              dat=data(sim_index).(var).Pxx(:,row);     
              legend_strings={'SUA','MUA'};
            case {'rastergram','raster'}
              set_name=regexp(var,'^([a-zA-Z0-9]+)_','tokens','once');
              allspikes{1}{1}=data(sim_index).([var '_spike_times']){row};
                % one pop, cell array of spike times for each cell in population
          end
          if num_rows>1
            text_string{row,col}=sprintf('cell %g',row);
          end
        % -----------------------------------------------------------------
        elseif num_sims==1 && num_pops==1 && num_vars>1
        % -----------------------------------------------------------------
          % one variable per row: dat = data(s=1).(var)(:,1:MTPP) where var=vars{v=r}
          var=var_fields{row};
          switch options.plot_type
            case 'waveform'
              dat=data(sim_index).(var);
            case 'power'
              AuxData=data(sim_index).([var '_Power_MUA']).Pxx;
              vlines=data(sim_index).([var '_Power_MUA']).PeakFreq;
              AuxDataName={'MUA Power'};
              var=[var '_Power_SUA'];
              dat=data(sim_index).(var).Pxx;
            case {'rastergram','raster'}
              set_name=regexp(var,'^([a-zA-Z0-9]+)_','tokens','once');
              allspikes{1}=data(sim_index).([var '_spike_times']);
          end
          shared_ylims_flag=0;
        % -----------------------------------------------------------------
        elseif num_sims==1 && num_pops>1 && num_vars==1
        % -----------------------------------------------------------------
          % one population per row: dat = data(s=1).(var)(:,1:MTPP) where var=vars{v=r}
          var=var_fields{row};
          switch options.plot_type
            case 'waveform'
              dat=data(sim_index).(var);
            case 'power'
              AuxData=data(sim_index).([var '_Power_MUA']).Pxx;
              vlines=data(sim_index).([var '_Power_MUA']).PeakFreq;
              AuxDataName={'MUA Power'};
              var=[var '_Power_SUA'];
              dat=data(sim_index).(var).Pxx;
            case {'rastergram','raster'}
              set_name=regexp(var,'^([a-zA-Z0-9]+)_','tokens','once');
              allspikes{1}=data(sim_index).([var '_spike_times']);
          end
        % -----------------------------------------------------------------
        elseif num_sims==1 && num_pops>1 && num_vars>1
        % -----------------------------------------------------------------
          % one population per row: dat = data(s=1).(var)(:,1:MTPP) where var=vars{these(p=r)}
          if isnan(pop_var_indices{row}(figset))
            continue;
          end
          var=var_fields{pop_var_indices{row}(figset)};
          switch options.plot_type
            case 'waveform'
              dat=data(sim_index).(var);         
            case 'power'
              AuxData=data(sim_index).([var '_Power_MUA']).Pxx;
              vlines=data(sim_index).([var '_Power_MUA']).PeakFreq;
              AuxDataName={'MUA Power'};
              var=[var '_Power_SUA'];
              dat=data(sim_index).(var).Pxx;
            case {'rastergram','raster'}
              set_name=regexp(var,'^([a-zA-Z0-9]+)_','tokens','once');
              allspikes{1}=data(sim_index).([var '_spike_times']);
          end
        % -----------------------------------------------------------------
        elseif num_sims>1 && num_pops==1 && num_vars==1
        % -----------------------------------------------------------------
          % one simulation per row: dat = data(s=r).(var)(:,1:MTPP) where var=vars{v=1}
          var=var_fields{1};
          switch options.plot_type
            case 'waveform'
              dat=data(sim_index).(var);
            case 'power'
              AuxData=data(sim_index).([var '_Power_MUA']).Pxx;
              vlines=data(sim_index).([var '_Power_MUA']).PeakFreq;
              AuxDataName={'MUA Power'};
              var=[var '_Power_SUA'];
              dat=data(sim_index).(var).Pxx;
            case {'rastergram','raster'}
              set_name=regexp(var,'^([a-zA-Z0-9]+)_','tokens','once');
              allspikes{1}=data(sim_index).([var '_spike_times']);
          end
        % -----------------------------------------------------------------
        elseif num_sims>1 && num_pops==1 && num_vars>1
        % -----------------------------------------------------------------
          % one simulation per row: dat = data(s=r).(var)(:,1:MTPP) where var=vars{v++}
          if isnan(pop_var_indices{1}(figset))
            continue;
          end
          var=var_fields{pop_var_indices{1}(figset)};
          switch options.plot_type
            case 'waveform'
              dat=data(sim_index).(var);              
            case 'power'
              AuxData=data(sim_index).([var '_Power_MUA']).Pxx;
              vlines=data(sim_index).([var '_Power_MUA']).PeakFreq;
              AuxDataName={'MUA Power'};
              var=[var '_Power_SUA'];
              dat=data(sim_index).(var).Pxx;
            case {'rastergram','raster'}
              set_name=regexp(var,'^([a-zA-Z0-9]+)_','tokens','once');
              allspikes{1}=data(sim_index).([var '_spike_times']);
          end
        % -----------------------------------------------------------------
        elseif num_sims>1 && num_pops>1 && num_vars==1
        % -----------------------------------------------------------------
          % one simulation per row, overlay pops: dat = <data(s=r).(var)(:,1:MTPP),2|vars>
          switch options.plot_type
            case 'waveform'
              % calculate averages across populations
              dat=nan(num_times,num_pops);
              for k=1:num_pops
                dat(:,k)=nanmean(data(sim_index).(var_fields{k}),2);
              end
              var=['<' variables{1} '>'];
            case 'power'
              dat=nan(length(xdata),num_pops);
              AuxData=nan(length(xdata),num_pops);
              AuxDataName={}; vlines=[];
              for k=1:num_pops
                dat(:,k)=nanmean(data(sim_index).([var_fields{k} '_Power_SUA']).Pxx,2);
                AuxData(:,k)=data(sim_index).([var_fields{k} '_Power_MUA']).Pxx;
                AuxDataName{end+1}=strrep([var_fields{k} '_Power_MUA'],'_','\_');
                vlines(end+1)=data(sim_index).([var_fields{k} '_Power_MUA']).PeakFreq;
              end
              var=['<' variables{1} '_Power_SUA>'];
            case {'rastergram','raster'}
              set_name={};
              for k=1:num_pops
                tmp=regexp(var_fields{k},'^([a-zA-Z0-9]+)_','tokens','once');
                set_name{k}=tmp{1};
                allspikes{k}=data(sim_index).([var_fields{k} '_spike_times']);
              end
              var=['<' variables{1} '>'];
              
          end
        % -----------------------------------------------------------------
        elseif num_sims>1 && num_pops>1 && num_vars>1
        % -----------------------------------------------------------------
          % one simulation per row, overlay pops: dat = <data(s=r).(var)(:,1:MTPP),2|vars(these)>
          switch options.plot_type
            case 'waveform'
              % calculate averages across populations
              dat=nan(num_times,num_pops);
              for k=1:num_pops
                if isnan(pop_var_indices{k}(figset))
                  continue;
                end
                var=var_fields{pop_var_indices{k}(figset)};
                dat(:,k)=nanmean(data(sim_index).(var),2);
              end
              var=['<' variables{figset} '>'];              
            case 'power'
              dat=nan(length(xdata),num_pops);
              AuxData=nan(length(xdata),num_pops);
              AuxDataName={}; vlines=[];
              for k=1:num_pops
                if isnan(pop_var_indices{k}(figset))
                  continue;
                end
                var=var_fields{pop_var_indices{k}(figset)};
                dat(:,k)=nanmean(data(sim_index).([var '_Power_SUA']).Pxx,2);
                AuxData(:,k)=data(sim_index).([var '_Power_MUA']).Pxx;
                AuxDataName{end+1}=strrep([var '_Power_MUA'],'_','\_');
                vlines(end+1)=data(sim_index).([var '_Power_MUA']).PeakFreq;
              end
              var=['<' variables{figset} '_Power_SUA>'];
            case {'rastergram','raster'}
              set_name={};
              for k=1:num_pops
                if isnan(pop_var_indices{k}(figset))
                  continue;
                end
                var=var_fields{pop_var_indices{k}(figset)};
                tmp=regexp(var,'^([a-zA-Z0-9]+)_','tokens','once');
                set_name{k}=tmp{1};
                allspikes{k}=data(sim_index).([var '_spike_times']);
              end
              var=['<' variables{figset} '>'];
          end
        end
        % #################################################################
        if size(dat,2)>1
          legend_strings=cellfun(@(x)['cell ' num2str(x)],num2cell(1:min(size(dat,2),max_legend_entries)),'uni',0);
        end
        
        if isfield(data,'varied')
          if num_sims>1
            % list the parameter varied along the rows first
            str=[row_param_name '=' num2str(row_param_values(row)) ': '];
            for k=1:num_varied
              fld=data(sim_index).varied{k};
              if ~strcmp(fld,row_param_name)
                val=data(sim_index).(fld);
                str=[str fld '=' num2str(val) ', '];
              end
            end
            if num_pops>1
              legend_strings=cellfun(@(x)[x ' (mean)'],pop_names,'uni',0);
            end
          else
            str='';
            for k=1:length(data.varied)
              fld=data(sim_index).varied{k};
              str=[str fld '=' num2str(data(sim_index).(fld)) ', '];
            end
          end
          text_string{row,col}=['(' strrep(str(1:end-2),'_','\_') ')'];
        end
        if ~isempty(AuxData) && length(legend_strings)<=max_legend_entries
          legend_strings=cat(2,legend_strings,AuxDataName);
        end
        % plot data
        %axes(haxes(axis_counter));
        set(gcf,'CurrentAxes',haxes(axis_counter));
        switch options.plot_type
          case {'waveform','power'}
            % finish preparing data
            if ~strcmp(options.yscale,'linear')
              dat=feval(options.yscale,dat); % log or log10
              % alternative approach: use semilogy for log10
            end
            if length(options.xlim)==2
              sel=(xdata>=options.xlim(1)&xdata<=options.xlim(2));
            else
              sel=1:length(xdata);
            end
            % plot traces
            if strcmp(options.plot_mode,'trace')
              % select max subset allowed
              dat=dat(:,1:min(size(dat,2),MTPP)); % select max subset to plot
              plot(xdata(sel),dat(sel,:));
            else
              imagesc(dat);
            end            
          case {'rastergram','raster'}
            % draw spikes
            ypos=0; % y-axis position tracker
            yticks=[]; % where to position population names
            yticklabels={}; % population names
            for p=1:length(allspikes) % loop over populations
              spikes=allspikes{p}; % spikes for one population
              for c=1:length(spikes) % loop over cells in population p
                spks=spikes{c}; % spikes for one cell
                for k=1:length(spks) % loop over spikes for cell c
                  spk=spks(k); % time of spike k in cell c of population p
                  line([spk spk],[c+ypos-.5 c+ypos+.5],'color','k'); hold on
                end
              end
              % record position for population tick name
              yticks(end+1)=ypos+c/2+.5; 
              yticklabels{end+1}=set_name{p};
              % draw line separating populations
              if length(allspikes)>1
                pos=c+ypos+.5;
                line([min(time) max(time)],[pos pos],'color','k','linewidth',3);
                if p<length(allspikes)
                  % increment y-position for next population
                  ypos=ypos+c;
                end
              end
            end
            % artificially set "dat" to get correct ylims below
            dat=[.5 ypos+c+.5];
            shared_ylims_flag=0;
            legend_strings='';
            % set y-ticks to population names
            set(gca,'ytick',yticks,'yticklabel',yticklabels);
            % set x-ticks
            plot([min(time) max(time)],[.5 .5],'w');
            nticks=length(get(gca,'xtick'));
            xticks=linspace(options.xlim(1),options.xlim(2),nticks);
            set(gca,'xtick',xticks,'xticklabel',xticks);
            %set(gca,'xticklabel',get(gca,'ytick'));
        end
        % plot auxiliary data
        if ~isempty(AuxData) %strcmp(options.plot_type,'power')
          hold on
          plot(xdata(sel),AuxData(sel,:),'-','linewidth',3);%,'o-','linewidth',3);%'--.');
        end
        % format axes
        if row==num_rows
          xlabel(xlab);
        else
          set(haxes(axis_counter),'XTickLabel',''); 
          %set(haxes(row),'YTickLabel','');
        end
        xlim(options.xlim);
        ylabel(strrep(var,'_','\_'));
        if ~isempty(options.ylim)
          ylims=options.ylim;
        elseif shared_ylims_flag
          % update max/min
          ylims(1)=min(ylims(1),min(dat(:)));
          ylims(2)=max(ylims(2),max(dat(:)));
        else
          % set ylim to max/min of this data set
          ylims=[min(dat(:)) max(dat(:))];
          if ylims(1)~=ylims(2)
            ylim(ylims);
          end
          % add text
          if ~isempty(text_string)
            xmin=min(xlim); xmax=max(xlim);
            ymin=min(dat(:)); ymax=max(dat(:));
            text_xpos=xmin+.05*(xmax-xmin);
            text_ypos=ymin+.9*(ymax-ymin);
            try
              text(text_xpos,text_ypos,text_string{row,col});
            end
          end  
        end
        % plot lines
        if ~isempty(vlines)
          for k=1:length(vlines)
            if ~isnan(vlines(k))
              line([vlines(k) vlines(k)],ylim,'color','k','linestyle','--');
            end
          end
        end
        % add legend
        if ~isempty(legend_strings) && axis_counter==1
          legend(legend_strings);
        end
      end % end loop over subplot columns
    end % end loop over subplot rows            
    % set y-limits to max/min over data in this figure
    if shared_ylims_flag || ~isempty(options.ylim)
      if ylims(1)~=ylims(2)
        set(haxes,'ylim',ylims);
      end
      if ~isempty(text_string)
        axis_counter=0;
        for row=1:num_rows
          for col=1:num_cols
            if ~ischar(text_string{row,col})
              continue;
            end
            axis_counter=axis_counter+1;
            %axes(haxes(axis_counter));
            set(gcf,'CurrentAxes',haxes(axis_counter));
            xmin=min(xlim); xmax=max(xlim);
            ymin=min(ylim); ymax=max(ylim);
            text_xpos=xmin+.05*(xmax-xmin);
            text_ypos=ymin+.9*(ymax-ymin);
            text(double(text_xpos),text_ypos,text_string{row,col});
          end
        end
      end
    end
  end % end loop over figures in this set
end % end loop over figure sets

% 1 sim, 1 pop, 1 var (X)
% 	N=1		one fig, one row (plot var X)
% 	N>1		one row per cell (plot var X), 			enough figs for all cells
% 	(nsims=1, num_pops=1, num_vars=1, pop_sizes>=1): num_fig_sets=1, num_figs=ceil(pop_sizes/MRPF), num_rows=min(pop_sizes,MRPF): row r: dat = data(s=1).(var)(:,c=r) where var=vars{v=1}
% 
% 1 sim, 1 pop, >1 vars (X,Y,...)
% 	N=1		one row per var (plot cell 1), 			enough figs for all vars
% 	N>1		one row per var (overlay cells), 		enough figs for all vars
% 	(nsims=1, num_pops=1, num_vars>=1, pop_sizes>=1): num_fig_sets=1, num_figs=ceil(num_vars/MRPF), num_rows=min(num_vars,MRPF): row r: dat = data(s=1).(var)(:,1:MTPP) where var=vars{v=r}
% 
% 1 sim, >1 pops, 1 var (X)
% 	all N=1		one row per pop (plot var X, cell 1),		enough figs for all pops
% 	any N>1		one row per pop (plot var X, overlay cells), 	enough figs for all pops
% 	(nsims=1, num_pops>=1, num_vars=1, pop_sizes>=1): num_fig_sets=1, num_figs=ceil(num_pops/MRPF), num_rows=min(num_pops,MRPF): row r: dat = data(s=1).(var)(:,1:MTPP) where var=vars{v=r}
% 
% 1 sim, >1 pops, >1 vars (X,Y,...)
% 	all N=1		one row per pop (plot var X, cell 1), 		enough figs for all pops, separate figs for each var
% 	any N>1		one row per pop (plot var X, overlay cells), 	enough figs for all pops, separate figs for each var
% 	(nsims=1, num_pops>=1, num_vars>=1, pop_sizes>=1): num_fig_sets=num_vars, num_figs=ceil(num_pops/MRPF), num_rows=min(num_pops,MRPF): row r: dat = data(s=1).(var)(:,1:MTPP) where var=vars{these(p=r)}
% 
% >1 sim, 1 pop, 1 var (X)
% 	N=1		one row per sim (plot var X, cell 1), 		enough figs for all sims
% 	N>1		one row per sim (plot var X, overlay cells), 	enough figs for all sims
% 	(nsims>1, num_pops=1, num_vars=1, pop_sizes>=1): num_fig_sets=1, num_figs=ceil(nsims/MRPF), num_rows=min(nsims,MRPF): row r: dat = data(s=r).(var)(:,1:MTPP) where var=vars{v=1}
% 
% >1 sim, 1 pop, >1 vars (X,Y,...)
% 	N=1		one row per sim (plot var X, cell 1), 		enough figs for all sims, separate figs for each var
% 	N>1		one row per sim (plot var X, overlay cells), 	enough figs for all sims, separate figs for each var
% 	(nsims>1, num_pops=1, num_vars=1, pop_sizes>=1): num_fig_sets=num_vars, num_figs=ceil(nsims/MRPF), num_rows=min(nsims,MRPF): row r: dat = data(s=r).(var)(:,1:MTPP) where var=vars{v++}
% 
% >1 sim, >1 pops, 1 var (X)
% 	all N=1		one row per sim (plot var X, overlay pops),	enough figs for all sims
% 	any N>1		one row per sim (plot var <X>, overlay pops),	enough figs for all sims
% 	(nsims>1, num_pops=1, num_vars=1, pop_sizes>=1): num_fig_sets=1, num_figs=ceil(nsims/MRPF), num_rows=min(nsims,MRPF): row r: dat = <data(s=r).(var)(:,1:MTPP),2|vars>
% 
% >1 sim, >1 pops, >1 vars (X,Y,...)
% 	all N=1		one row per sim (plot var X, overlay pops),	enough figs for all sims, separate figs for each var
% 	any N>1		one row per sim (plot var <X>, overlay pops),	enough figs for all sims, separate figs for each var
% 	(nsims>1, num_pops=1, num_vars=1, pop_sizes>=1): num_fig_sets=num_vars, num_figs=ceil(nsims/MRPF), num_rows=min(nsims,MRPF): row r: dat = <data(s=r).(var)(:,1:MTPP),2|vars(these)>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
