function model=CombineModels(model1,model2)
%% model=CombineModels(model1,model2)
% Purpose: combine subfields in two DynaSim model structures
% Inputs: two models to be combined
% Output: DynaSim model with fields combined from both models
% 
% See also: CheckModel, GenerateModel

% standardize model structures
model1=CheckModel(model1);
model2=CheckModel(model2);

% combine fields from sub-structures
model.parameters=concatenate_structures(model1.parameters,model2.parameters);
model.fixed_variables=concatenate_structures(model1.fixed_variables,model2.fixed_variables);
model.functions=concatenate_structures(model1.functions,model2.functions);
model.monitors=concatenate_structures(model1.monitors,model2.monitors);
model.ODEs=concatenate_structures(model1.ODEs,model2.ODEs);
model.ICs=concatenate_structures(model1.ICs,model2.ICs);

% concatenate cell and structure arrays
model.state_variables=cat(2,model1.state_variables,model2.state_variables);
model.conditionals=cat(2,model1.conditionals,model2.conditionals);
model.linkers=cat(2,model1.linkers,model2.linkers);
model.comments=cat(2,model1.comments,model2.comments);

% combine .specification from model1 and model2 (this is necessary for
% building a new model from two indepedent models to which connection
% mechanisms are added...)
% todo: call something like old combine_models() function from old DynaSim
% ...

% standardize resulting model
% model=CheckModel(model);
  % note: if this call to CheckModel() is uncommented-out, the changes noted
  % in CheckModel() should also be made...

% reorder fields according to first input
model=orderfields(model,model1);

% SUBFUNCTIONS
function out=concatenate_structures(a,b)
if isempty(a) && ~isempty(b)
  out=b;
elseif ~isempty(a) && isempty(b)
  out=a;
elseif isempty(a) && isempty(b)
  out=a;
elseif ~isempty(a) && ~isempty(b)
  out=catstruct(a,b);
end
