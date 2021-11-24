function eqns=dsExtractModelStrings(MODEL,display_mode,display_flag)
% Purpose: construct string to display DynaSim model equations:
%   ODEs
%   ICs
%   conditionals
%   functions
%   fixed_variables
%   parameters
%   monitors
%   comments
%
% display_mode {'model' (default),'specification','odefun','xpp'}
% 
% Example: display DynaSim model equations
% dsExtractModelStrings(MODEL,'model',1);
%
% Example: display and use script defining DynaSim specification structure
% eqns=dsExtractModelStrings(MODEL,'specification',1);
% spec=eval(eqns); % note: equivalent to MODEL.specification
%
% Example: integrate system using built-in Matlab solver
%   eqns=dsExtractModelStrings(MODEL,'odefun',0);
%   fun=eval(eqns{1});
%   ic=eqns{2};
%   [t,y]=ode23(fun,[0 100],ic);
%   figure; plot(t,y);

% display_mode OPTIONS (todo: add options 2 and 3):
% 1. Display resulting model equations from DynaSim model structure
% 2. Display ODEFUN (function handle string for ODE system: @(X,t)...)
%     tip: use fun=str2func(eqns{1}) to obtain function handle from output
%          and ic=eval(eqns{2}) to obtain initial condition vector.
% 3. Display script defining the DynaSim specification structure
% 4. Display XPP .ode model implementation (see notes below)

if nargin<3
  display_flag=0;
end
if nargin<2 || isempty(display_mode)
  display_mode='model';
end

eqns={};
switch lower(display_mode)
  case 'model' % Display resulting model equations from DynaSim model structure
    % standardize DynaSim model structure
    MODEL=dsCheckModel(MODEL);
    % ODEs and ICs:
    if ~isempty(MODEL.state_variables)
      eqns{end+1}='% DIFFERENTIAL EQUATIONS:';
      vars=MODEL.state_variables;
      for i=1:length(vars)
        eqns{end+1}=sprintf('%%  %s'' = %s',vars{i},MODEL.ODEs.(vars{i}));
      end
      eqns{end+1}='%';
      eqns{end+1}='% Initial conditions:';
      for i=1:length(vars)
        eqns{end+1}=sprintf('%%  %s(0) = %s',vars{i},MODEL.ICs.(vars{i}));
      end
      eqns{end+1}='';
    end
    % conditionals
    if ~isempty(MODEL.conditionals)
      eqns{end+1}='% CONDITIONALS:';
      for i=1:length(MODEL.conditionals)
        str=sprintf('  if(%s)then(%s)',MODEL.conditionals(i).condition,MODEL.conditionals(i).action);
        if ~isempty(MODEL.conditionals(i).else)
          str=sprintf('%selse(%s)',str,MODEL.conditionals(i).else);
        end
        eqns{end+1}=sprintf('\t%s',str);
      end
      eqns{end+1}='';
    end
    types={'parameters','fixed_variables','functions'};%,'monitors'
    type_headers={'% PARAMETERS:','% FIXED VARIABLES:','% FUNCTIONS:','% MONITORS:'};
    for p=1:length(types)
      type=types{p};
      header=type_headers{p};
      if ~isempty(MODEL.(type))
        eqns{end+1}=header;
        fields=fieldnames(MODEL.(type));
        for i=1:length(fields)
          val=MODEL.(type).(fields{i});
          if ~ischar(val)
            val=toString(val,'compact');
          end
          eqns{end+1}=sprintf('%s = %s',fields{i},val);
        end
      end
      eqns{end+1}='';
    end
  case 'odefun' % Display ODEFUN (function handle string for ODE system: @(X,t)...)
    % Approach:
    % 1. evaluate params -> fixed_vars -> funcs
    % 2. evaluate ICs to get (# elems) per state var
    % 3. prepare state vector X
    % 4. replace state vars in ODEs by X
    % 5. combine X ODEs into ODEFUN

    % evaluate params -> fixed_vars -> funcs
    types={'parameters','fixed_variables','functions'};
    for p=1:length(types)
      type=types{p};
      if ~isempty(MODEL.(type))
        fields=fieldnames(MODEL.(type));
        for i=1:length(fields)
          val=MODEL.(type).(fields{i});
          if ~ischar(val)
            val=toString(val,'compact');
          end
          % evaluate
          eval(sprintf('%s = %s;',fields{i},val));
        end
      end
    end

    % evaluate ICs to get (# elems) per state var and set up generic state var X
    num_vars=length(MODEL.state_variables);
    num_elems=zeros(1,num_vars);
    old_vars=MODEL.state_variables;
    new_vars=cell(1,num_vars);
    new_inds=cell(1,num_vars);
    all_ICs=cell(1,num_vars);
    IC_names={};
    state_var_index=0;
    for i=1:num_vars
      var=MODEL.state_variables{i};
      % evaluate ICs to get (# elems) per state var
      ic=eval([MODEL.ICs.(var) ';']);
      num_elems(i)=length(ic);
      % set state var indices a variables for generic state vector X
      all_ICs{i}=ic;
      IC_names{i}=repmat({var},[1 num_elems(i)]);
      new_inds{i}=state_var_index+(1:length(ic));
      new_vars{i}=sprintf('X(%g:%g)',new_inds{i}(1),new_inds{i}(end));
      state_var_index=state_var_index+length(ic);
    end

    % prepare ODE system (comma-separated ODEs)
    ODEs=strtrim(struct2cell(MODEL.ODEs));
    idx=cellfun(@isempty,regexp(ODEs,';$')); % lines that need semicolons
    ODEs(idx)=cellfun(@(x)[x ';'],ODEs(idx),'uni',0);
    ODEs=[ODEs{:}]; % concatenate ODEs into a single string
    ODEs=strrep(ODEs,';',','); % replace semicolons by commas

    % substitute in generic state vector X
    for i=1:num_vars
      ODEs=dynasim_strrep(ODEs,old_vars{i},new_vars{i});
    end

    % prepare outputs (function handle string, ICs, and element names for
    % mapping each X(i) to a particular state variable):
    ODEFUN = eval(['@(t,X) [' ODEs '];']);
    IC=cat(2,all_ICs{:});
    elem_names=cat(2,IC_names{:});

    eqns{1}=ODEFUN;
    eqns{2}=IC;
    eqns{3}=elem_names;

    %{
      % usage:

      eqns=dsExtractModelStrings(MODEL,'odefun',0);
      ODEFUN=eqns{1};
      IC=eqns{2};
      elem_names=eqns{3};

      dt=.01; t=0:dt:100;
      y=zeros(length(t),length(IC));
      y(1,:)=IC;
      for i=2:length(t)
        y(i,:)=y(i-1,:)+dt*ODEFUN(t,y(i-1,:));
      end
      figure; plot(t,y); legend(elem_names{:},'Location','EastOutside');

      y=IC;
      for i=1:1e4
        y=y+dt*ODEFUN(0,y);
      end;

    %}

  otherwise
    error('options ''specification'' and ''xpp'' not implemented yet.');
end

if display_flag
  cellfun(@disp,eqns);
end
