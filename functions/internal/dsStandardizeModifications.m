function modifications = dsStandardizeModifications(modifications,specification, varargin)
% dsStandardizeModifications - convert all modifications into 3-column cell
%                              matrix format (namespace,variable,value) and expand
%
% dsStandardizeModifications expands multiple namespace in one cell to multiple
% cells, ensures source->target direction, moves specific mech to variable name from
% namespace, and converts mech dot notation to underscores
%
% Inputs:
%   modifications: modifications structure
%
% Inputs (optional):
%   specification: needed if namespace empty to check for first population name.
%                  just need "specification.populations(1).name)".

% Dev Notes:
%   E.R. Feb 2018: Moved to separate function

if nargin < 2
  specification = [];
end

%% auto_gen_test_data_flag argin
options = dsCheckOptions(varargin,{'auto_gen_test_data_flag',0,{0,1}},false);
if options.auto_gen_test_data_flag
  varargs = varargin;
  varargs{find(strcmp(varargs, 'auto_gen_test_data_flag'))+1} = 0;
  varargs(end+1:end+2) = {'unit_test_flag',1};
  argin = [{modifications},{specification}, varargs]; % specific to this function
end

if isstruct(modifications)
  % TODO
  % convert structure to cell matrix
  % ...
end

% check backward compatibility
modifications=backward_compatibility(modifications);

% check for empty object specification
missing_objects=find(cellfun(@isempty,modifications(:,1)));
if any(missing_objects)
  if isempty(specification)
    error('Need non-empty specification input for population name')
  end
  
  % set to default population name
  for i=1:length(missing_objects)
    modifications{missing_objects(i),1}=specification.populations(1).name;%'pop1';
  end
end

% support modifying multiple elements (of namespace or variable) simultaneously
% approach -- add extra entry for each thing to modify
% eg) expand {'(E,I)','Cm',2} to {'E','Cm',2; 'I','Cm',2}
% note: should be able to support {'(E,I)','(EK,EK2)',-80}
if any(~cellfun(@isempty,regexp(modifications(:,1),'^\(.*\)$'))) || ...
    any(~cellfun(@isempty,regexp(modifications(:,2),'^\(.*\)$')))
  
  % loop over modifications
  modifications_={};
  for i=1:size(modifications,1)
    % check namespace for ()
    namespaces=regexp(modifications{i,1},'[\w\.\-<>]+','match');
    
    % check variable for ()
    variables=regexp(modifications{i,2},'[\w\.-]+','match');
    
    if ischar(modifications{i,3})
      
      % expand list of modifications
      for j=1:length(namespaces)
        for k=1:length(variables)
          modifications_(end+1,1:3)={namespaces{j},variables{k},modifications{i,3}};
        end
      end
      
    elseif isnumeric(modifications{i,3})
      
      % check size of values matches number of namespaces, variables
      if isscalar(modifications{i,3}) % in case number of values is one
        modifications{i,3} = repmat(modifications{i,3},length(variables),length(namespaces));
      else
        if size(modifications{i,3},1) ~= length(variables) || size(modifications{i,3},2) ~= length(namespaces)
          % in case values is number of variables x 1
          if size(modifications{i,3},1) == length(variables) && size(modifications{i,3},2) == 1
            modifications{i,3} = repmat(modifications{i,3},1,length(namespaces));
            % in case values is 1 x number of namespaces
          elseif size(modifications{i,3},2) == length(namespaces) && size(modifications{i,3},1) == 1
            modifications{i,3} = repmat(modifications{i,3},length(variables),1);
            % TODO: char inputs
            % elseif ischar(modifications{i,3})
            % string input
          else % if ~ischar(modifications{i,3})
            error(['Numerical values varied over must be in array format,',...
              'where dimensions 1, 2, and 3 correspond to mechanisms, values, and populations varied over.'])
          end
        end
      end
      
      % expand list of modifications
      for j=1:length(namespaces)
        for k=1:length(variables)
          modifications_(end+1,1:3)={namespaces{j},variables{k},modifications{i,3}(k,j)};
        end
      end
      
    end % ischar
    
  end % mods
  
  modifications=modifications_;
  
end % if groups

% standardize connection object
for i = 1:size(modifications,1)
  obj=modifications{i,1}; % population name or connection source-target
  fld=modifications{i,2}; % population name, size, or parameter name
  
  % convert target<-source to source->target
  if any(strfind(obj,'<-'))
    if any(obj=='.') % check for dot reference to connection mechanism
      ind=find(obj=='.');
      suffix=obj(ind:length(obj));
      obj=obj(1:ind-1);
    else
      suffix='';
    end
    ind=strfind(obj,'<-');
    obj=[obj(ind(1)+2:end) '->' obj(1:ind(1)-1)];
    obj=[obj suffix];
    
    % update modifications
    modifications{i,1} = obj; % population name or connection source-target
  end
  
  % check for dot notation. Change POP.MECH to MECH_PARAM notation.
  if any(obj=='.')
    tmp=regexp(obj,'\.','split');
    obj=tmp{1};
    MECH=tmp{2};
    fld=[MECH '_' fld];
    
    % update modifications
    modifications{i,1} = obj; % population name or connection source-target
    modifications{i,2} = fld; % population name, size, or parameter name
  end
  
  % check for dot notation. Change MECH.PARAM to MECH_PARAM notation.
  if any(fld=='.')
    tmp=regexp(fld,'\.','split');
    fld=[tmp{1} '_' tmp{2}];
    
    % update modifications
    modifications{i,2} = fld; % population name, size, or parameter name
  end
end

%% auto_gen_test_data_flag argout
if options.auto_gen_test_data_flag
  argout = {modifications}; % specific to this function
  
  dsUnitSaveAutoGenTestDataLocalFn(argin, argout); % localfn
end

end


%% Local Fn

function modifications = backward_compatibility(modifications)
% convert 2-column specification to 3-column specification with empty object name
if size(modifications,2)==2
  tmp={};
  for i=1:size(modifications,1)
    tmp{i,1}='';
    tmp{i,2}=modifications{i,1};
    tmp{i,3}=modifications{i,2};
  end
  modifications=tmp;
end

% convert 4-column specification to 3-column
if size(modifications,2)==4
  for i=1:size(modifications,1)
    if strcmp(modifications{i,2},'parameters')
      % shift parameter name to 2nd column
      modifications{i,2}=modifications{i,3};

      % shift parameter value to 3rd column
      modifications{i,3}=modifications{i,4};
    end
  end

  % remove fourth column
  modifications=modifications(:,1:3);
end

% convert connection reference source-target to source->target
if any(~cellfun(@isempty,regexp(modifications(:,1),'\w-\w')))
  % find elements to adjust
  inds=find(~cellfun(@isempty,regexp(modifications(:,1),'\w-\w')));
  for i=1:length(inds)
    modifications{inds(i),1}=strrep(modifications{inds(i),1},'-','->');
  end
end

end