function data = probeFI(model,varargin)
%% data=ds.probeFI(model,varargin)
% purpose: run experiment delivering varying levels of tonic input to each 
% population across multiple simulations.
% inputs:
%   model - DynaSim model structure
%   options:
%     target - ...
%     amplitudes - ...
%     (any options for ds.simulateModel)
%     (any options for ds.calcFR)
% outputs:
%   array of DynaSim data structures for simulations varying inputs
% 
% see also: ds.simulateModel, ds.calcFR

options=ds.checkOptions(varargin,{...
  'target','ODE1',[],...
  'amplitudes',0:2:10,[],...
  'onset',0,[],...
  },false);

model=ds.checkModel(model);

npops=length(model.specification.populations);
modifications=cell(npops,3);
vary=cell(npops,3);
for i=1:npops
  name=model.specification.populations(i).name;
  % prepare list of modifications to add tonic drive to all populations in model
  %modifications(i,:)={name,'equations',['cat(' options.target ',+TONIC)']};
  modifications(i,:)={name,'equations',sprintf('cat(%s,+TONIC*(t>%g))',options.target,options.onset)};
  % prepare specification to vary tonic drive across simulations
  vary(i,:)={name,'TONIC',options.amplitudes};
end

% apply modifications
model=ds.applyModifications(model,modifications);

% run simulations varying tonic drives
data=ds.simulateModel(model,'vary',vary,varargin{:});

% prepare ds.calcFR options:
if ~any(cellfun(@(x)isequal(x,'variable'),varargin))
  % default to extract spike events from the first state variable
  var=regexp(model.state_variables{1},'_(\w+)$','tokens','once');
  if ~isempty(var)
    varargin{end+1}='variable';
    varargin{end+1}=['*_' var{1}];
  end
end
% calculate resulting firing rates for each population
data=ds.calcFR(data,varargin{:});
