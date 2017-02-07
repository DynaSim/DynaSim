function [data,studyinfo] = ImportData(file,varargin)
%% [data,studyinfo]=ImportData(data_file)
% data=ImportData(data_file)
% Purpose: load data into DynaSim formatted data structure.
% Inputs:
%   data_file - data file name in accepted format (csv, mat, ...), or
%     also accepted: list of data files, studyinfo structure, study_dir, or studyinfo file
%   studyinfo - DynaSim studyinfo structure (see CheckStudyinfo)
%   options -
%     'verbose_flag',1,{0,1},...
%     'process_id',[],[],... % process identifier for loading studyinfo if necessary
%     'time_limits',[],[],... % [beg,end] ms (see NOTE 2)
%     'variables',[],[],... % cell array of matrix names (see NOTE 2)
% Outputs:
%   DynaSim data structure:
%     data.labels           : list of state variables and monitors recorded
%     data.(state_variables): state variable data matrix [time x cells]
%     data.(monitors)       : monitor data matrix [time x cells]
%     data.time             : time vector [time x 1]
%     data.simulator_options: simulator options used to generate simulated data
%     data.model            : model used to generate simulated data
%     [data.varied]         : list of varied model components
%     [data.results]        : list of derived data sets created by post-processing
% 
% NOTE 1: CSV file structure
%   assumes CSV file contains data organized according to output from
%   WriteDynaSimSolver: time points along rows; state variables and 
%   monitors are columns; first column is time vector; next columns are
%   state variables; final columns are monitors. first row has headers
%   for each column. if a population has more than one cell, different
%   cells are sequential columns with same header repeated for each cell.
% 
% NOTE 2: DynaSim data exported to MAT-files are HDF-compatible. To obtain
% partial data sets without having to load the entire file, use ImportData
% with options 'time_limits' and/or 'variables'. Alternatively, the entire
% data set can be loaded using ImportData with default options, then
% subsets extracted using SelectData with appropriate options.
% 
% Examples:
% 1) full data set
% data=ImportData('data.mat'); % load single data set
% data=ImportData(studyinfo); % load all data sets in studyinfo.study_dir
% 2) partial data set with HDF-style loading
% data=ImportData('data.mat','variables','pop1_v','time_limits',[1000 4000])
% 
% See also: SimulateModel, ExportData, CheckData, SelectData

% todo: specify subsets to return in terms of varied parameters, time_limits, ROIs, etc
% possible format for specifying range_varied:
% {'E','gNa',[.1 .3]; 'I->E','tauI',[15 25]; 'I','mechanism_list','+iM'}
% idea: achieve by calling function SelectData() at end of this function.

% Check inputs
options=CheckOptions(varargin,{...
  'verbose_flag',1,{0,1},...
  'process_id',[],[],... % process identifier for loading studyinfo if necessary
  'time_limits',[],[],...
  'variables',[],[],...
  },false);

if ischar(options.variables)
  options.variables={options.variables};
end
% check if input is a DynaSim studyinfo structure
if ischar(file) && isdir(file) % study directory
  study_dir=file;
  clear file
  file.study_dir=study_dir;
end
if isstruct(file) && isfield(file,'study_dir')
  % "file" is a studyinfo structure.
  % retrieve most up-to-date studyinfo structure from studyinfo.mat file
  studyinfo=CheckStudyinfo(file.study_dir,'process_id',options.process_id);
  % get list of data_files from studyinfo
  data_files={studyinfo.simulations.data_file};
  success=cellfun(@exist,data_files)==2;
  data_files=data_files(success);
  sim_info=studyinfo.simulations(success);
  % load each data set recursively
  keyvals=Options2Keyval(options);
  num_files=length(data_files);
  for i=1:num_files
    fprintf('loading file %g/%g: %s\n',i,num_files,data_files{i});
    tmp_data=ImportData(data_files{i},keyvals{:});
    num_sets_per_file=length(tmp_data);
    if ~isfield(tmp_data,'varied')
    % add varied info 
      % this is necessary here when loading .csv data lacking metadata
      tmp_data.varied={};
      modifications=sim_info(i).modifications;
      modifications(:,1:2) = cellfun( @(x) strrep(x,'->','_'),modifications(:,1:2),'UniformOutput',0);
      for j=1:size(modifications,1)
        varied=[modifications{j,1} '_' modifications{j,2}];
        for k=1:num_sets_per_file
          tmp_data(k).varied{end+1}=varied;
          tmp_data(k).(varied)=modifications{j,3};
        end
      end    
    end
    % store this data
    if i==1
      total_num_sets=num_sets_per_file*num_files;
      set_indices=0:num_sets_per_file:total_num_sets-1;
      % preallocate full data matrix based on first data file
      data(1:total_num_sets)=tmp_data(1);
%       data(1:length(data_files))=tmp_data;
%     else
%       data(i)=tmp_data;
    end
    % replace i-th set of data sets by these data sets
    data(set_indices(i)+(1:num_sets_per_file))=tmp_data;
  end
  return;
else
  studyinfo=[];
end
% check if input is a list of data files (todo: eliminate duplicate code by
% combining with the above recursive loading for studyinfo data_files)
if iscellstr(file)
  data_files=file;
  success=cellfun(@exist,data_files)==2;
  data_files=data_files(success);
  keyvals=Options2Keyval(options);
  % load each data set recursively
  for i=1:length(data_files)
    tmp_data=ImportData(data_files{i},keyvals{:});
    % store this data
    if i==1
      % preallocate full data matrix based on first data file
      data(1:length(data_files))=tmp_data;
    else
      % replace i-th data element by this data set
      data(i)=tmp_data;
    end
  end
  return;  
end

if ischar(file)
  [~,~,ext]=fileparts(file);
  switch lower(ext)
    case '.mat'
      % MAT-file contains data fields as separate variables (-v7.3 for HDF)
      if isempty(options.time_limits) && isempty(options.variables)
        % load full data set
        data=load(file);
        % if file only contains a structure called 'data' then return that
        if isfield(data,'data') && length(fieldnames(data))==1
          data=data.data;
        end
      else
        % load partial data set
        % use matfile() to load HDF subsets given varargin options...
        obj=matfile(file); % MAT-file object
        varlist=who(obj); % variables stored in mat-file
        labels=obj.labels; % list of state variables and monitors
        if iscellstr(options.variables) % restrict variables to load
          labels=labels(ismember(labels,options.variables));
        end
        simulator_options=obj.simulator_options;
        time=(simulator_options.tspan(1):simulator_options.dt:simulator_options.tspan(2))';
        time=time(1:simulator_options.downsample_factor:length(time));
        if ~isempty(options.time_limits)
          % determine time indices to load
          time_indices=nearest(time,options.time_limits(1)):nearest(time,options.time_limits(2));
        else
          % load all time points
          time_indices=1:length(time);
        end
        % create DynaSim data structure:
        data=[];
        data.labels=labels;
        % load state variables and monitors
        for i=1:length(labels)
          data.(labels{i})=obj.(labels{i})(time_indices,:);
        end
        data.time=time(time_indices);
        data.simulator_options=simulator_options;
        if ismember('model',varlist)
          data.model=obj.model;
        end
        if ismember('varied',varlist)
          varied=obj.varied;
          data.varied=varied;
          for i=1:length(varied)
            data.(varied{i})=obj.(varied{i});
          end
        end
        if ismember('results',varlist)
          results=obj.results;
          if iscellstr(options.variables)
            results=results(ismember(results,options.variables));
          end
          data.results=results;
          % load results
          for i=1:length(results)
            data.(results{i})=obj.(results{i})(time_indices,:);
          end
        end
      end
    case '.csv'
      % assumes CSV file contains data organized according to output from WriteDynaSimSolver:
      data=ImportCSV(file);
      if ~(isempty(options.time_limits) && isempty(options.variables))
        % limit to select subsets
        data=SelectData(data,varargin{:}); % todo: create SelectData()
      end
    otherwise
      error('file type not recognized. ImportData currently supports DynaSim data structure in MAT file, data values in CSV file.');
  end
end
