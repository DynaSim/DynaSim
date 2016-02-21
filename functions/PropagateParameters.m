function model=PropagateParameters(model,varargin)
%% model = SubstituteParameters(model)
% purpose: substitute parameter values or prepend parameter names with 
%          prefix across all model equations.
% input: DynaSim model structure
% options:
%   action: {'substitute','prepend'} (default: substitute)
%   prefix: string prepended to all parameter names if action is 'prepend'
% output: DynaSim model structure with updated parameter in all equations
% 
% see also: PropagateFunctions, WriteDynaSimSolver

% Check inputs
model=CheckModel(model);
if ~isstruct(model.parameters)
  % nothing to do
  return;
end
% Check inputs
options=CheckOptions(varargin,{...
  'action','substitute',{'substitute','prepend'},...
  'prefix','pset.p.',[],...
  },false);

parameters=model.parameters;

%% 1.0 Propagate through sub-structures
target_types={'fixed_variables','functions','monitors','ODEs','ICs'};
% loop over types of model data
for type_index=1:length(target_types)
  type=target_types{type_index};
  % info for this type
  s=model.(type);
  if isstruct(s)
    update_these=fieldnames(s);
    expressions=struct2cell(s);
    % loop over target expressions from which to eliminate internal function calls
    for i=1:length(expressions)    
      if isempty(expressions{i})
        continue;
      end
      % update expressions of this type
      switch options.action
        case 'substitute'
          expressions{i}=insert_parameters(expressions{i},parameters,[]);
        case 'prepend'
          expressions{i}=insert_parameters(expressions{i},parameters,options.prefix);
      end
    end
    % update model with expressions that have parameter values in them
    model.(type)=cell2struct(expressions,update_these,1);
  end
end

%% 2.0 Propagate parameters through structure arrays (conditionals)
if ~isempty(model.conditionals)
  target_types={'condition','action','else'};
  for type_index=1:length(target_types)
    type=target_types{type_index};
    expressions={model.conditionals.(type)};
    % loop over conditional expressions from which to eliminate internal function calls
    for i=1:length(expressions)
      if isempty(expressions{i})
        continue; 
      end
      % update expressions of this type
      switch options.action
        case 'substitute'
          expressions{i}=insert_parameters(expressions{i},parameters,[]);
        case 'prepend'
          expressions{i}=insert_parameters(expressions{i},parameters,options.prefix);
       end     
    end
    [model.conditionals(1:length(model.conditionals)).(type)]=deal(expressions{:});
  end
end

function expression=insert_parameters(expression,parameters,prefix)
  if isnumeric(expression)
    % convert to string and return string
    expression=toString(expression);
    return;
  end
  allwords=regexp(expression,'[a-zA-Z]+\w*','match');
  words=unique(allwords);
  found_parameters=words(ismember(words,fieldnames(parameters)));
  if ~isempty(found_parameters)
    % substitute those found into this target expression
    for ff=1:length(found_parameters)
      % name of found parameter
      found_parameter=found_parameters{ff};
      if isempty(prefix) % no prefix given, substitute value instead
        % found value to replace found parameter name in target
        found_value=parameters.(found_parameter);
        % convert found value into string
        if isnumeric(found_value)
          if length(found_value)>1
            found_value=sprintf('[%s]',num2str(found_value));
          else
            found_value=num2str(found_value);
          end
        elseif iscell(found_value)
          if iscellstr(found_value)
            tmp=cellfun(@(x)['''' x '''' ','] ,found_value,'uni',0);
          else
            tmp=cellfun(@(x)[num2str(x) ','] ,found_value,'uni',0);
          end
          tmp=[tmp{:}];
          found_value=sprintf('{%s}',tmp(1:end-1));
        elseif isa(found_value,'function_handle')
          found_value=func2str(found_value);
        end
      else % prefix provided, substitute prefix_name
        found_value=[prefix found_parameter];
      end
      if ~ischar(found_value)
        warning('failed to convert parameter ''%s'' to string and substitute into model equations:',found_parameter);
        found_value
      else
        % update expression
        num_found = length(find(ismember(allwords,found_parameter)));
        for iter=1:num_found
          expression=dynasim_strrep(expression,found_parameter,found_value);
        end
      end
    end
  end

