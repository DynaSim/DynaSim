function spec=CheckSpecification(specification)
%% specification=CheckSpecification(specification)
% Purpose: standardize specification structure and auto-populate missing fields
% Input: DynaSim specification structure or equations
% Output: DynaSim specification structure (standardized)
% 
% DynaSim specification structure
% .populations(i) (required): contains info for defining independent
%     population models
%   .name (default: 'pop1')      : name of population
%   .size (default: 1)           : number of elements in population (i.e., # cells)
%   .equations (required)        : string listing equations (see NOTE 1)
%   .mechanism_list (default: []): cell array listing mechanisms (see NOTE 2)
%   .parameters (default: [])    : parameters to assign across all equations in
%     the population. provide as cell array list of key/value pairs
%     {'param1',value1,'param2',value2,...}
%   .conditionals (default: [])  : (see NOTE 3) if-then conditional actions
%     .namespace (auto)    : string giving the namespace of the condition 
%                            (pop_ or pop_mech_)
%     .condition (required): string giving the condition to check
%     .action (required)   : what to do if the condition is met
%     .else (default: [])  : what to do if the condition is not met
%   .monitors (default: [])      : (see NOTE 3) substructure with fields specifying what to
%     record on each step of numerical integration in addition to state variables.
%   .model (default: [])   : optional DynaSim model structure
% .connections(i) (default: []): contains info for linking population models
%   .source (required if >1 pops): name of source population (see NOTE 7)
%   .target (required if >1 pops): name of target population
%   .mechanism_list (required)   : list of mechanisms that link two populations
%   .parameters (default: [])    : parameters to assign across all equations in
%     mechanisms in this connection's mechanism_list.
% 
% NOTE 1: .equations can be an equation string, cell array listing
% equation strings, or a file name pointing to a model / equations stored
% on disk (accepted file types: .eqns (equations of population), 
% .m (function defining a model structure), ...)
% 
% NOTE 2: .mechanism_list is a cell array listing names of mechanisms to be
% included in the population or used to connect two populations. each
% mechanism name must have a mechanism file with the same name somewhere in
% the search path (the file should have extension .mech). The search path
% starts with the current directory, then the subdirectories of
% [dynasim]/models, and lastly the full matlab search path.
% 
% NOTE 3: conditionals and monitors are most easily specified by including
% them in .equations.
% e.g.) spec.populations.equations='dv/dt=-v; if(v<eps)(v=10); monitor o=v^2'
%       data=SimulateModel(spec); figure; plot(data.time,data.pop1_o);
% 
% NOTE 4: "pops" can be used instead of "populations". "cons" can be used
% instead of "connections".
% 
% NOTE 5: all population info can be embedded in the equation string.
% Specify name by starting the string with 'NAME: *' (e.g., 'E: dv/dt=-v').
% Specify size by including [SIZE] after the state variable (e.g., 'dv[5]/dt=-v').
% Specify mechanism_list by including cell array listing mechanism names
% without single quotes (e.g., 'dv/dt=@current; {iNa,iK}').
% 
% NOTE 6: the mechanism linker target IDENTIFIER used to link mechanism variables and
% functions to population equations can be overriden by including
% @NEWIDENTIFIER in the equations and after the mechanism name. (e.g.,
% 'dv/dt=@M; {iNa,iK}@M'; or .mechanism_list={'iNa@M','iK@M'}).
% 
% NOTE 7: "direction" can be used instead of "source" and "target". The
% syntax is "SOURCE->TARGET" or "TARGET<-SOURCE". Either will be properly
% split into .source and .target fields. SOURCE and TARGET must be existing
% population names.
% 
% Example 1: obtain empty specification structure with all fields
% specification=CheckSpecification([]);
% 
% Example 2: standardize existing specification
% specification=CheckSpecification(specification)
% 
% Example 3: standardize equations in cell array
% eqns={
%   's=10; r=27; b=2.666';
%   'dx/dt=s*(y-x)';
%   'dy/dt=r*x-y-x*z';
%   'dz/dt=-b*z+x*y';
% specification=CheckSpecification(eqns);
% 
% Example 4: standardize equations in character array
% eqns='tau=10; R=10; E=-70; dV/dt=(E-V+R*1.55)/tau; if(V>-55)(V=-75)';
% specification=CheckSpecification(eqns);
% 
% Example 5: standardize specification with compact field names
% s.pops.size=10;
% s.pops.equations='dv/dt=-v';
% s.cons.mechanism_list='iGABAa';
% s=CheckSpecification(s)
% 
% Example 6: standardize specification with everything in equation string
% s.pops.equations='E:dv[10]/dt=@M+I; {iNa,iK}@M; I=10';
% s=CheckSpecification(s)
% 
% see also: GenerateModel, CheckModel

% check if input is a string or cell with equations and package in spec structure
if ischar(specification) || iscell(specification)
  spec.populations.equations=specification;
elseif isstruct(specification)
  spec=specification;
elseif isempty(specification)
  spec=struct;
else
  error('specification must be a DynaSim specification structure or a string with equations or sub-model filename.');
end

spec=backward_compatibility(spec);
pop_field_order={'name','size','equations','mechanism_list','parameters',...
  'conditionals','monitors','model'};
con_field_order={'source','target','mechanism_list','parameters'};

if ~isfield(spec,'populations')
  spec.populations.name='pop1';
end
if ~isfield(spec,'connections')
  spec.connections=[];
end

% 1.0 standardize populations
if ~isfield(spec.populations,'name')
  spec.populations(1).name='pop1';
end
if ~isfield(spec.populations,'size')
  spec.populations(1).size=1;
end
if ~isfield(spec.populations,'equations')
  spec.populations(1).equations=[];
end
if ~isfield(spec.populations,'mechanism_list')
  spec.populations(1).mechanism_list=[];
end
if ~isfield(spec.populations,'parameters')
  spec.populations(1).parameters={};
end
if ~isfield(spec.populations,'conditionals')
  spec.populations(1).conditionals=[];
end
if ~isfield(spec.populations,'monitors')
  spec.populations(1).monitors=[];
end
if ~isfield(spec.populations,'model')
  spec.populations(1).model=[];
end
% special case: split equations with '[...][...]...[...]' into multiple populations
for i=1:length(spec.populations)
  eqn=spec.populations(i).equations;
  if ~isempty(eqn) && ischar(eqn)
    if ~isempty(regexp(eqn,'\[[a-z_A-Z].*\]','match','once'))
      % create extra population
      tmp=regexp(eqn(2:end-1),'\],?\s*\[','split');
      spec.populations(i).equations=tmp{1};
      for j=2:length(tmp)
        spec.populations(end+1)=spec.populations(i);
        spec.populations(end).equations=tmp{j};
        spec.populations(end).name=sprintf('pop%g',length(spec.populations));
      end
    end
  end
end
% standardize each population separately
for i=1:length(spec.populations)
  % population names
  if isempty(spec.populations(i).name)
    spec.populations(i).name=sprintf('pop%g',i);
  end
  % population sizes
  if isempty(spec.populations(i).size)
    spec.populations(i).size=1;
  end
  % make mechanism list a cell array of mechanism names
  if ischar(spec.populations(i).mechanism_list)
    spec.populations(i).mechanism_list={spec.populations(i).mechanism_list};
  end
  % parameter cell arrays
  if ~iscell(spec.populations(i).parameters)
    spec.populations(i).parameters={};
  end
  % standardize equations
  if ~isempty(spec.populations(i).equations)
    % convert cell array of equations into character array
    if iscell(spec.populations(i).equations)
      eqns=spec.populations(i).equations;
      for k=1:length(eqns)
        eqn=eqns{k};
        if ~isempty(eqn) && ~strcmp(eqn(end),';')
          eqns{k}(end+1)=';';
        end
      end
      spec.populations(i).equations=[eqns{:}];
    end
    % extract name from equations (e.g., TC:...)
    eqn=spec.populations(i).equations;
    name=regexp(eqn,'^\w+:','match','once');
    if ~isempty(name)
      % remove name indicator from equation
      eqn=strrep(eqn,name,'');
      % store name in specification
      name=regexp(name,'^(\w+):','tokens','once');
      spec.populations(i).name=name{1};
      spec.populations(i).equations=eqn;
    end
    % extract size from equations if present (eg, v[4]'=.., dv[4]/dt=...)
    eqn=spec.populations(i).equations;
    pattern='((\w+(\[\d+\])'')|(d\w+(\[\d+\])/dt))\s*='; % support size spec, dv[4]/dt
    % extract all differentials with size specification
    LHSs=regexp(eqn,pattern,'match');
    if ~isempty(LHSs)
      % extract sizes from all differentials (eg, 4 from v[4] or dv[4]/dt)
      szs=nan(1,length(LHSs));
      for k=1:length(LHSs)
        tmp=regexp(LHSs{k},'\w+\[(\d+)\]''','tokens','once');
        if isempty(tmp)
          tmp=regexp(LHSs{k},'d\w+\[(\d+)\]/dt','tokens','once');
        end
        szs(k)=str2num(tmp{1});
        % remove size from ODE in population equations
        old=LHSs{k};
        new=strrep(LHSs{k},['[' tmp{1} ']'],'');
        eqn=strrep(eqn,old,new);
      end
      % check that all vars in same population have same size
      if ~all(szs==szs(1))
        error('all variables in same population must have same size. split ODEs with different sizes into different populations.');
      end
      spec.populations(i).equations=eqn;
      spec.populations(i).size=szs(1);
    end
    % add mechanisms embedded in equations to mechanism_list ({M1,M2,...})
    % ----------
    % todo: make the following a subfunction and apply it also to connection
    % mechanism lists (eg, for cons.mechanism_list='{AMPA,NMDA}@M')
    % ----------
    % extract mechanism list from equations
    %mech_lists=regexp(spec.populations(i).equations,'\s*{.*}\s*;?\s*','match');
    %mech_lists=regexp(spec.populations(i).equations,'\s*(\w+:)?{.*}\s*(@\w+)?;?\s*','match');
     mech_lists=regexp(spec.populations(i).equations,'\s*(\w+:)?{[\w\d@:,]*}\s*(@\w+)?;?\s*','match');
    % test: mech_list=regexp('v''=@M+sin(2*pi*t); {iNa, iK}','{.*}','match');
    if ~isempty(mech_lists)
      for k=1:length(mech_lists)
        mech_list=strtrim(mech_lists{k});
        % remove mechanism list from equations
        spec.populations(i).equations=strtrim(strrep(spec.populations(i).equations,mech_list,''));
        % append external link alias to each internal mechanism name (eg, {a,b}@M, alias @M)
        external_link=regexp(mech_list,'}(@[\w\d]+;?)','tokens');        
        if ~isempty(external_link)
          % get external link alias
          external_link=[external_link{:}];
          % remove external link alias from mech_list
          mech_list=strrep(mech_list,external_link{1},'');
          % remove ';' from alias before appending to mech names
          external_link=strrep(external_link{1},';',''); 
          % get list of mechanism names in cell array
          %words=regexp(mech_list,'[a-zA-Z]+[\w:]*','match');
          words=regexp(mech_list(2:end-1),',','split');
          % append external link alias to each mechanism name
          for w=1:length(words)
            mech_list=strrep(mech_list,words{w},[words{w} external_link]);
          end
        end
        % prepend host name to each internal mechanism name (eg,
        % infbrain:{a,b} -> {infbrain:a,infbrain:b}
        host_name=regexp(mech_list,';?\s*([\w\d]+):{','tokens','once'); 
        if ~isempty(host_name)
          % get external link alias
          host_name=[host_name{:}];
          % remove external link alias from mech_list
          mech_list=strrep(mech_list,[host_name ':'],'');
          % get list of mechanism names in cell array
          words=regexp(mech_list(2:end-1),',','split');
          % append external link alias to each mechanism name
          for w=1:length(words)
            mech_list=strrep(mech_list,words{w},[host_name ':' words{w}]);
          end
        end
        % split into list of mechanism names
        mechanisms=regexp(mech_list,'[\w:@]+','match');        
        % append mechanism from equations to mechanism_list
        if iscell(spec.populations(i).mechanism_list)
          spec.populations(i).mechanism_list=cat(2,mechanisms,spec.populations(i).mechanism_list);
        else
          spec.populations(i).mechanism_list=mechanisms;
        end
      end
    end
    % extract population-level parameters from equations
    eqn=spec.populations(i).equations;
    p=getfield(ParseModelEquations(eqn),'parameters');
    if ~isempty(p)
      param_name=fieldnames(p);
      param_value=struct2cell(p);
      for l=1:length(param_name)
        value=eval(param_value{l});
        if isempty(spec.populations(i).parameters)
          spec.populations(i).parameters={param_name{l},value};
        elseif ~ismember(param_name{l},spec.populations(i).parameters(1:2:end))
          spec.populations(i).parameters{end+1}=param_name{l};
          spec.populations(i).parameters{end+1}=value;
        end
      end
    end
    % incorporate user-supplied parameters in pop equations if used in them
    if ~isempty(spec.populations(i).parameters)
      keys=spec.populations(i).parameters(1:2:end);
      vals=spec.populations(i).parameters(2:2:end);
      % add user-supplied params to pop equations if present in them
      % approach: look for populations.parameters in population.equations that are not explicitly 
      % defined in population.equations and append their definition explicitly to pop.eqns
      eqn=spec.populations(i).equations;
      % get list of parameters/variables/functions in population equations
      words=unique(regexp(eqn,'[a-zA-Z]+\w*','match'));
      % find those in user-supplied parameters
      found_words=words(ismember(words,keys));
      if ~isempty(found_words)     
        % set in population equations if not already defined there
        for ff=1:length(found_words)
          found_word=found_words{ff};
          % check if not explicitly set in population equations
          if isempty(regexp(eqn,[';\s*' found_word '\s*='],'once')) && ... % not in middle or at end
             isempty(regexp(eqn,['^' found_word '\s*='],'once')) % not at beginning
            % explicitly set in population equations
            if eqn(end)~=';', eqn(end+1)=';'; end % add semicolon if necessary
            precision=8; % number of digits allowed for user-supplied values
            found_value = toString(vals{strcmp(found_word,keys)},precision);
            eqn=[eqn sprintf(' %s=%s;',found_word,found_value)];
          end
        end
        spec.populations(i).equations=eqn;
      end
    end
  end
  % expand mechanism list if any element is itself a list of mechanisms (eg, {'iCa','{CaBuffer,iCan}'} or '{CaBuffer,iCan}')
  spec.populations(i).mechanism_list=expand_list(spec.populations(i).mechanism_list);  
end

% 2.0 standardize connections
if ~isempty(spec.connections)
  % check for proper fields in connections substructure
  if ~isfield(spec.connections,'source')
    spec.connections(1).source=[];
  end
  if ~isfield(spec.connections,'target')
    spec.connections(1).target=[];
  end
  if ~isfield(spec.connections,'mechanism_list')
    spec.connections(1).mechanism_list=[];
  end
  if ~isfield(spec.connections,'parameters')
    spec.connections(1).parameters={};
  end  
end
for i=1:length(spec.connections)
  if isempty(spec.connections(i).source) && length(spec.populations)==1
    spec.connections(i).source=spec.populations(1).name;
    spec.connections(i).target=spec.populations(1).name;
  elseif isempty(spec.connections(i).source) && length(spec.populations)>1
    error('connection source and target populations must be specified in specification.connections when the model contains more than one population.');
  end
  % make mechanism list a cell array of mechanism names
  if ischar(spec.connections(i).mechanism_list)
    spec.connections(i).mechanism_list={spec.connections(i).mechanism_list};
  end
  % expand mechanism list if any element is itself a list of mechanisms (eg, {'AMPA','{GABAa,GABAb}'} or '{GABAa,GABAb}')
  spec.connections(i).mechanism_list=expand_list(spec.connections(i).mechanism_list);
  % parameter cell arrays
  if ~iscell(spec.connections(i).parameters)
    spec.connections(i).parameters={};
  end
end

% remove populations with size==0
sizes=[spec.populations.size];
if any(sizes==0)
  % find null populations
  null_pops=find(sizes==0);
  null_names={spec.populations(null_pops).name};
  % remove from .populations
  spec.populations(null_pops)=[];
  % remove from connections
  if ~isempty(spec.connections)
    sources={spec.connections.source};
    targets={spec.connections.target};
    null_conns=ismember(sources,null_names) | ismember(targets,null_names);
    spec.connections(null_conns)=[];
  end
end    

% 3.0 sort fields
% remove extra fields
otherfields=setdiff(fieldnames(spec.populations),pop_field_order);
spec.populations=rmfield(spec.populations,otherfields);
% sort standardized fields
spec.populations=orderfields(spec.populations,pop_field_order);
if isstruct(spec.connections)
  otherfields=setdiff(fieldnames(spec.connections),con_field_order);
  spec.connections=rmfield(spec.connections,otherfields);
  spec.connections=orderfields(spec.connections,con_field_order);
end

% 4.0 replace mechanism names by full file names
% this is necessary so that regenerated models will use the same mechanism
% files to recreate the model (e.g., when a cluster job simulates a
% modified version of an original base model).
[~,files]=LocateModelFiles(spec);
if ~isempty(files)
  fnames={};
  for f=1:length(files)
    [~,name]=fileparts(files{f});
    fnames{f}=name;
  end
  % update population and connection mechanism lists
  fields={'populations','connections'};
  for f=1:length(fields)
    object=fields{f};
    for i=1:length(spec.(object))
      for j=1:length(spec.(object)(i).mechanism_list)
        mech=spec.(object)(i).mechanism_list{j};
        if ismember(mech,fnames)
          spec.(object)(i).mechanism_list{j}=files{find(ismember(fnames,mech),1,'first')};
        end
      end
    end
  end
end

function list=expand_list(list)
% expand mechanism list if any element is itself a list of mechanisms (eg, {'AMPA','{GABAa,GABAb}'} or '{GABAa,GABAb}')
if isempty(list)
  return;
end
if any(~cellfun(@isempty,regexp(list,'[{,}]+')))
  mechs={};
  for k=1:length(list)
    tmp=regexp(list{k},'\w+','match');
    mechs=cat(2,mechs,tmp{:});
  end
  list=mechs;
end

function spec=backward_compatibility(spec)
% purpose: change name of fields from old to new convention
% rename "nodes" or "entities" to "populations"
if isfield(spec,'nodes')
  spec.populations=spec.nodes;
  spec=rmfield(spec,'nodes');
end
if isfield(spec,'cells')
  spec.populations=spec.cells;
  spec=rmfield(spec,'cells');
end
if isfield(spec,'entities')
  spec.populations=spec.entities;
  spec=rmfield(spec,'entities');
end
if isfield(spec,'pops')
  spec.populations=spec.pops;
  spec=rmfield(spec,'pops');
end
if isfield(spec,'cons')
  spec.connections=spec.cons;
  spec=rmfield(spec,'cons');
end
if isfield(spec,'populations')
  % rename population "label" to "name"
  if isfield(spec.populations,'label')
    for i=1:length(spec.populations)
      spec.populations(i).name=spec.populations(i).label;
    end
    spec.populations=rmfield(spec.populations,'label');
  end
  % rename population "multiplicity" to "size"
  if isfield(spec.populations,'multiplicity')
    for i=1:length(spec.populations)
      spec.populations(i).size=spec.populations(i).multiplicity;
    end
    spec.populations=rmfield(spec.populations,'multiplicity');
  end
  % rename population "dynamics" to "equations"
  if isfield(spec.populations,'dynamics')
    for i=1:length(spec.populations)
      spec.populations(i).equations=spec.populations(i).dynamics;
    end
    spec.populations=rmfield(spec.populations,'dynamics');
  end
  % rename population "mechanisms" to "mechanism_list"
  if isfield(spec.populations,'mechanisms')
    for i=1:length(spec.populations)
      spec.populations(i).mechanism_list=spec.populations(i).mechanisms;
    end
    spec.populations=rmfield(spec.populations,'mechanisms');
  end
end
% check for old (pre,post) organization of connections substructure
if isfield(spec,'connections') && size(spec.connections,1)>1
  % convert to linear connections structure array
  old=spec.connections;
  spec=rmfield(spec,'connections');
  index=1;
  for i=1:size(old,1)
    for j=1:size(old,2)
      if ~isempty(old(i,j).mechanisms)
        spec.connections(index).source=spec.populations(i).name;
        spec.connections(index).target=spec.populations(j).name;
        if isfield(old,'mechanisms')
          spec.connections(index).mechanism_list=old(i,j).mechanisms;
        elseif isfield(old,'mechanism_list')
          spec.connections(index).mechanism_list=old(i,j).mechanism_list;
        end
        if isfield(old,'parameters')
          spec.connections(index).parameters=old(i,j).parameters;
        end
        index=index+1;
      end
    end
  end
end
if isfield(spec,'connections') && isfield(spec.connections,'direction')
  for i=1:length(spec.connections)
    if ischar(spec.connections(i).direction)
      str=spec.connections(i).direction;
      if any(regexp(str,'->','once'))
        pops=regexp(str,'->','split');
        spec.connections(i).source=pops{1};
        spec.connections(i).target=pops{2};
      elseif any(regexp(str,'<-','once'))
        pops=regexp(str,'<-','split');
        spec.connections(i).source=pops{2};
        spec.connections(i).target=pops{1};
      end
    end
  end
  spec.connections=rmfield(spec.connections,'direction');
end