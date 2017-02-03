function [paths,files]=LocateModelFiles(input)
%% [paths,files]=LocateModelFiles(input)
% Purpose: locate mechanism files associated with DynaSim specifications.
% Input: DynaSim specification or model structure or string or cell array of
%        strings listing mechanism names or files.
% Outputs: 
%   paths - unique paths to mechanism files
%   files - full names of files containing mechanism sub-models
% 
% See also (used by): ParseModelEquations, CheckHostPaths, CreateBatch

% extract list of mechanisms from input
if ischar(input)
  % convert to cell array of strings
  mechanism_list={input};
elseif iscellstr(input)
  mechanism_list=input;
elseif isstruct(input)
  if isfield(input,'specification')
    % extract specification from DynaSim model structure
    input=input.specification;
  elseif isfield(input,'model') && isfield(input.model,'specification') 
    % extract specification from DynaSim data structure
    input=input.model.specification;
  elseif isfield(input,'base_model') && isfield(input.base_model,'specification') 
    % extract specification from DynaSim studyinfo structure
    input=input.base_model.specification;
  else
    % this is probably a specification structure
  end
  mechanism_list={};
  if isfield(input,'populations') && isstruct(input.populations)
    % extract mechanism_list from populations in DynaSim specification structure
    m={input.populations.mechanism_list};
    m=unique([m{:}],'stable');
    mechanism_list=cat(2,mechanism_list,m);
    % add equations to mechanism_list in case it contains a .eqns file
    for i=1:length(input.populations)
      if ~isempty(input.populations(i).equations)
        mechanism_list{end+1}=input.populations(i).equations;
      end
    end
  end
  if isfield(input,'connections') && isstruct(input.connections)
    % extract mechanism_list connections in DynaSim specification structure
    m={input.connections.mechanism_list};
    m=unique([m{:}],'stable');
    mechanism_list=cat(2,mechanism_list,m);
  end  
end

% remove @ pointers from mechanism identifiers
if any(~cellfun(@isempty,regexp(mechanism_list,'@','once')))
  mechanism_list=regexp(mechanism_list,'^([^@]+)@?','tokens','once');
  mechanism_list=[mechanism_list{:}];
end
% exclude elements with non-word characters (these are not file names)
keep=cellfun(@isempty,regexp(mechanism_list,'[^\w\.\-/]'));
%keep=cellfun(@isempty,regexp(mechanism_list,'[^\w\.]'));
mechanism_list=mechanism_list(keep);

% locate mechanism files
files={};
if iscellstr(mechanism_list)
  for f=1:length(mechanism_list) % loop over mechanisms
    % determine name of mechanism file (assuming recommended extensions)
    mech=mechanism_list{f};
    if exist(mech,'file')
      file=mech;
    elseif exist([mech '.eqns'],'file')
      file=[mech '.eqns'];
    elseif exist([mech '.mech'],'file')
      file=[mech '.mech'];
    elseif exist([mech '.txt'],'file')
      file=[mech '.txt'];
    elseif exist([mech '.m'],'file')
      file=[mech '.m'];
    else
      % TODO: Add error catching here!
      file='';
    end
    % use 'which' to get full filename of file in Matlab path
    if isempty(fileparts(file)) % file is a name without a path
      file=which(file);
    end
    if ~isempty(file)
      files{end+1}=file;
    end
  end
end

% extract unique paths to mechanism files
paths=unique(cellfun(@fileparts,files,'uni',0),'stable');
