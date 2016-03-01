function [model,name_map] = ParseModelEquations(text,varargin)
%% model = ParseModelEquations(STRING,'param',value,...)
% Purpose: parse equations and organize model data in DynaSim model structure
% Inputs: 
%   STRING (required): one of --
%   - string with equations
%   - string with name of file containing equations (.eqns or .mech)
%   options (using key/value pairs: 'option1',value1,...):
%   - 'namespace': (default: '') added as prefix to beginning of parameter/etc names
%   - 'delimiter': (default: ';') separates expressions on same line of model text
%   user-supplied parameter values: ('key',value): name (key) of parameters to be set and associated user-supplied values
% Outputs:
%   model: DynaSim model structure (see CheckModel for details)
%   name_map: {name, namespace_name, namespace, type}, useful for namespace-specific substitutions across multiple sub-models
%             (see description in GenerateModel for more information)
% 
% NOTE 1: .eqns files contain fully self contained model equations; 
% .mech files define (sub)models that depend on variables linked from 
% elsewhere. However, this function does not distinguish between the two.
% 
% Examples:
% model = ParseModelEquations('dx/dt=3*a*x; x(0)=0','a',0);
% model = ParseModelEquations('dx/dt=3*a*x, x(0)=0','a',0,'delimiter',',');
% model = ParseModelEquations('CalciumPump.mech','namespace','HH');
% model = ParseModelEquations('LIFneuron.eqns');
% model = ParseModelEquations('a=2; b=2*a; f(x)=b; dx/dt=f(x); x(0)=0; if(x>1)(x=0); current->f(x); monitor f(x); % comments')
% 
% parsing individual sub-models from specification:
% equations=specification.populations(1).equations;
% [model,map] = ParseModelEquations(equations,'namespace','pop')
% population_mechanism=specification.populations(1).mechanism_list{1};
% [model,map] = ParseModelEquations(population_mechanism,'namespace','pop_mech')
% connection_mechanism=specification.connections(1).mechanism_list{1};
% [model,map] = ParseModelEquations(connection_mechanism,'namespace','pop_pop_mech')
% 
% See also: ClassifyEquation, GenerateModel, LocateModelFiles

model=[];
name_map={};

% organize optional user-supplied info
% key/value pairs
if nargin>2 % check for at least text input and one key/value pair
  keys=varargin(1:2:end); % parameters to set
  values=varargin(2:2:end); % values to use
else
  keys=[];
  values=[];
end
% set namespace
if ismember('namespace',keys) % check for user-supplied namespace (i.e., namespace)
  namespace=values{ismember(keys,'namespace')}; % user-supplied namespace
  if ~isempty(namespace)
    namespace=[namespace '_'];
  else
    namespace='';
  end
else
  namespace='';
end
% set delimiter
if ismember('delimiter',keys) % check for user-supplied delimiter
  delimiter = values(ismember(keys,'delimiter')); % user-supplied delimiter
else
  delimiter=';';
end
% error handling for improper input format
if ~ischar(namespace)
  error('model "namespace" must be a string.');
end
if ~ischar(delimiter)
  error('expression "delimiter" must be a string.');
end

%% 1.0 convert text into cell array of strings (one string per line)
% check for DynaSim extensions if input is single \w+ string
if ischar(text) && ~any(which(text)) && isempty(regexp(text,'[^\w.]','once')) % isempty(regexp(text,'[^\w]','once'))
%if ischar(text) && ~exist(text,'file') && isempty(regexp(text,'[^\w.]','once')) % isempty(regexp(text,'[^\w]','once'))
  [~,text]=LocateModelFiles(text);
  if iscell(text) && ~isempty(text)
    text=text{1};
  end
end
% check if input is a filename
if ischar(text) && exist(text,'file')
  [~,name,ext]=fileparts(text);
  switch ext
    case '.m'
      model=feval(name); % evaluate model-creating function and return model
      return;
    case '.mat' % todo: uncomment once ImportModel supports loading .mat
      %model=ImportModel(text);
      %return;
  end
  % load equations from file
  [text,res]=readtext(text,'\n','%'); % text: cell array of strings, one element per line in text file
  % remove all lines without text
  text=text(res.stringMask);
  % remove leading/trailing white space
  text=strtrim(text);
  % end each line with semicolon
  for i=1:length(text)
    if ~isequal(text{i}(end),';')
      text{i}(end+1)=';';
    end
  end
  % concatenate into a single string
  text=[text{:}]; % concatenate text from all lines  
end

% split string into cell array of lines delimited by semicolon
if ischar(text)
  % remove end-line semicolon if present so split lines are free of all ';'
  if text(end)==';'
    text=text(1:end-1);
  end
  % account for the one exception where ';' does not delimit lines: 
  % conditional actions with multiple statements (expr1; expr2)
  % approach: replace ';' by ',' here then reverse the replacement below
  % when storing the action in model.conditionals
  pattern='(if\([^;]+\)\s*\([^;\)]+);([^;]+\))'; % if(condiiton)(action1;action2)
  replace='$1,$2';
  text=regexprep(text,pattern,replace,'ignorecase');
  % now split string into cell array of lines
  text = strtrim(regexp(text,delimiter,'split'));
end
if ~iscellstr(text)
  error('input not recognized. equations must be provided in single string, cell array of strings, or a text file');
end

%% 2.0 classify and parse lines; store info in model structure
model.parameters=[];
model.fixed_variables=[];
model.functions=[];
model.monitors=[];
model.state_variables={};
model.ODEs=[];
model.ICs=[];
model.conditionals=[];
model.linkers=[];
model.comments={};
for index=1:length(text) % loop over lines of text
  % organize model data in model structure
  line=text{index}; % choose a single expression (lines with multiple expressions have already been split above using regexp-split)
  [line,comment]=remove_comment(line); % remove comments
  if isempty(line) % e.g., entire line was a comment, or there was nothing there originally
    if ~isempty(comment)
      model.comments{end+1}=comment;
    end
    continue;
  end
  switch ClassifyEquation(line,delimiter) % classify
    case 'parameter'        % var=(string or number)
      rhs=regexp(line,'=(.+)$','tokens','once');
      lhs=regexp(line,'^([\w\.]+)\s*=','tokens','once');
      lhs{1}=strrep(lhs{1},'.','_'); % e.g., Na.g --> Na_g
      name=strtrim(lhs{1}); expression=rhs{1};
      model.parameters.([namespace name]) = expression;
      name_map(end+1,:) = {name,[namespace name],namespace,'parameters'};
      if ~isempty(comment)
        model.comments{end+1}=sprintf('%s (parameter): %s',[namespace name],comment);
      end
    case 'fixed_variable'   % var=(expression with grouping or arithmetic), var(#), var([#]), var([# #]), var([#,#]), var(#:#), var(#:end), var([#:#]), var([#:end])
      lhs=regexp(line,'^(\w+)\s*=','tokens','once');
      rhs=regexp(line,'=(.+)$','tokens','once');
      name=strtrim(lhs{1}); expression=rhs{1};
      model.fixed_variables.([namespace name]) = expression;
      name_map(end+1,:) = {name,[namespace name],namespace,'fixed_variables'};
      if ~isempty(comment)
        model.comments{end+1}=sprintf('%s (fixed_variable): %s',[namespace name],comment);
      end
    case 'function'         % f(vars)=exression
%       if any(line=='@')
%         line=strrep(line,'@','');
%         %error('model specification error: delete the ''@'' character from all function definitions and try again.');
%       end
      name=regexp(line,'^(.+)\(.*\)\s*=','tokens','once');
      vars=regexp(line,'\((.+)\)\s*=','tokens','once');      
      rhs=regexp(line,'=(.+)$','tokens','once');
      name=strtrim(name{1}); 
      expression=sprintf('@(%s)%s',vars{1},rhs{1});
      model.functions.([namespace name]) = expression;
      name_map(end+1,:) = {name,[namespace name],namespace,'functions'};
      if ~isempty(comment)
        model.comments{end+1}=sprintf('%s (function): %s',[namespace name],comment);
      end
    case 'ODE'              % x'=expression or dx/dt=expression
      var=regexp(line,'^d(\w+)/dt\s*=','tokens','once'); % x from dx/dt=
      if isempty(var)
        var=regexp(line,'^(\w+)''\s*=','tokens','once'); % x from x'=
      end
      rhs=regexp(line,'=(.+)$','tokens','once');
      state_variable=strtrim(var{1}); expression=rhs{1};
      model.ODEs.([namespace state_variable])=expression;      
      if ~ismember([namespace state_variable],model.state_variables)
        name_map(end+1,:) = {state_variable,[namespace state_variable],namespace,'state_variables'};
        model.state_variables{end+1}=[namespace state_variable];
      end
      if ~isempty(comment)
        model.comments{end+1}=sprintf('d/dt %s (ODE): %s',[namespace state_variable],comment);
      end
    case 'IC'               % x(0)=expression
      var=regexp(line,'^(\w+)\(','tokens','once');
      rhs=regexp(line,'=(.+)$','tokens','once');
      state_variable=strtrim(var{1}); expression=rhs{1};
      model.ICs.([namespace state_variable])=expression;
      if ~isempty(comment)
        model.comments{end+1}=sprintf('%s(0) (IC): %s',[namespace state_variable],comment);
      end
    case 'monitor'          % monitor f=(expression or function)
      % split list of monitors
      lines=strtrim(regexp(line,',','split')); 
      % loop over monitors in list
      for l=1:length(lines)
        % process this monitor
        line=lines{l};
        % split left and right parts of monitor
        lhs=regexp(line,'^monitor ([\w,@\s\.]+)','tokens','once');
        if isempty(lhs)
          lhs=regexp(line,'([\w,@\s\.]+)','tokens','once');
        end        
        rhs=regexp(line,'=(.+)$','tokens','once');      
        % expand list of monitor names (e.g., monitor iNa.I, iK.I)
        names=strtrim(regexp(lhs{1},',','split'));
        for i=1:length(names) % loop over list of monitors on this line
          name=names{i};
          % process special monitors (those including '.', e.g., v.spikes(0))
          % todo: clean up or generalize this procedure...
          if any(name=='.')
            % check for numeric monitor argument
            arg=regexp(line,[name '\(([-+]*\w+)\)'],'tokens','once');
            %arg=regexp(line,[name '\(([-+]*\d+)\)'],'tokens','once');
            % set argument as expression (see WriteDynaSimSolver() for usage as such)
            if ~isempty(arg)
              rhs=arg;
            end
          end
          % convert into valid monitor name
          name=strrep(name,'.','_'); % index sub-namespace (monitor Na.I)        
          if ~isempty(rhs), expression=rhs{1}; else expression=[]; end
          model.monitors.([namespace name]) = expression;
          name_map(end+1,:) = {name,[namespace name],namespace,'monitors'};
          if ~isempty(comment)
            model.comments{end+1}=sprintf('%s (monitor): %s',[namespace name],comment);
          end
        end
      
      end
      
    case 'conditional'      % if(conditions)(actions)
      groups=regexp(line,')(','split');
      condition=regexp(groups{1},'^if\s*\((.*)','tokens','once');
      if length(groups)==2
        if groups{2}(end)==')'
          groups{2}=groups{2}(1:end-1);
        end
        then_action=groups{2};
        else_action=[];
      elseif numel(groups==3)
        if groups{3}(end)==')'
          groups{3}=groups{3}(1:end-1);
        end
        then_action=groups{2};
        else_action=groups{3};
      end
      model.conditionals(end+1).namespace=namespace;
      model.conditionals(end).condition=condition{1};
      model.conditionals(end).action=strrep(then_action,',',';'); % restore semicolon-delimited multiple actions like if(x>1)(x=0;y=0)
      if length(groups)>2
        model.conditionals(end).else=else_action;
      else
        model.conditionals(end).else=[];
      end
%       groups=regexp(line,'\(((\w+\([\w,@]+\))?[\w@-\+*^/\s><=,&|\.]+)\)','tokens');
%       condition=groups{1}{1};
%       model.conditionals(end+1).namespace=namespace;
%       model.conditionals(end).condition=condition;
%       model.conditionals(end).action=strrep(groups{2}{1},',',';'); % restore semicolon-delimited multiple actions like if(x>1)(x=0;y=0)
%       if length(groups)>2
%         model.conditionals(end).else=groups{3}{1};
%       else
%         model.conditionals(end).else=[];
%       end
      if ~isempty(comment)
        model.comments{end+1}=sprintf('%s conditional(%s): %s',namespace,condition,comment);
      end
    case 'linker'           % [link ]? target operation expression (e.g., link target += f(x))
      % viable options: ((\+=)|(-=)|(\*=)|(/=)|(=>))
      line=regexprep(line,'^link (\s*\w)','$1');
      lhs=regexp(line,'^([^\+\-*/=]+)','tokens','once'); % +=
      rhs=regexp(line,'=>?(.+)$','tokens','once');
      model.linkers(end+1).namespace=namespace;
      if ~isempty(lhs), target=strtrim(lhs{1}); else target=[]; end
      if ~isempty(rhs), expression=strtrim(rhs{1}); else expression=[]; end
      if expression(end)==';', expression=expression(1:end-1); end
      if isempty(target), target=expression; end % for sharing state var across mechanisms in same population
      if isempty(expression), expression=target; end
      model.linkers(end).target=target;
      model.linkers(end).expression=expression;
      model.linkers(end).operation='+=';
      if ~isempty(comment)
        model.comments{end+1}=sprintf('%s linkers(%s->%s): %s',namespace,target,expression,comment);
      end
      if ~isempty(comment)
        model.linkers(end).comment=comment;
      end
    otherwise
      warning('ignoring line, failed to classify :%s',line);
  end
end


% subfunctions
function [line,comment]=remove_comment(line)
% purpose: split line into model content and comment
index=find(line=='%',1,'first');
if isempty(index) 
  % check other valid comment delimiters
  index=find(line=='#',1,'first');
end
if isempty(index) % no comment found
  comment='';
else % split line into comment and non-comment line sections
  comment=line(index:end);
  line=line(1:index-1);
end
