function data = imdpLoadDSdata(data_dir, varargin)
if ~exist('data_dir', 'var') || isempty(data_dir)
  data_dir = pwd;
end

% if ~exist('options', 'var') || isempty(options)
%   options = struct();
% end

%% Check Options
options = CheckOptions(varargin,{...
  'classifyFn', [], [],...
  'overwriteBool', 0, {0,1},...
  'verboseBool', 1, {0,1},...
  },false);


%% Load or Create imdpData
filePath = fullfile(data_dir, 'imdpData.mat');
if ~exist(filePath,'file') || options.overwriteBool
  % Determine classifyFn
  if isempty(options.classifyFn)
    dataList = lscell(fullfile(data_dir, 'data'));
    classFnStrList = regexpi(dataList, '_(class.*)\.mat', 'tokens');
    classFnStrList = classFnStrList(~cellfun(@isempty,classFnStrList));
    classFnStrList = cellfun(@(x) x{1}, classFnStrList);
    classFnStrList = unique(classFnStrList);
    if numel(classFnStrList) == 1
      classFnStr = classFnStrList{1};
      classifyFn = str2func(classFnStr);
    else
      error('Found multiple class functions. Specifiy the desired one in options.classifyFn')
    end
  end
  
  % Import DS Data
  vfprintf('Importing SaveVaried results...\n')
  VariedData = ImportResults(data_dir, @SaveVaried);
  vfprintf('\tDone importing SaveVaried results.\n')
  
  vfprintf('Importing classification results...\n')
  classResults = ImportResults(data_dir, classifyFn);
  vfprintf('\tDone importing classification results\n')
  
  if length(classResults) ~= length(VariedData)
    wprintf('Different lengths for SaveVaried results and %s results. Results may not be accurate.', char(classifyFn))
  end
  
  vfprintf('Preparing data to save...\n')
  
  simIDs = num2cell((1:length(VariedData))', 2);
  
  % Remove missing data
  missingVariedDataInd = any(cellfun(@isempty,struct2cell(VariedData)));
  missingClassResultsInd = cellfun(@isempty,classResults);
  VariedData(missingVariedDataInd) = [];
  classResults(missingVariedDataInd) = [];
  simIDs(missingVariedDataInd) = [];
  VariedData(missingClassResultsInd) = [];
  classResults(missingClassResultsInd) = [];
  simIDs(missingClassResultsInd) = [];
  
  % Get class info
  try
    info = feval(classifyFn, 'info');
    Classes = struct();
    Classes.names = info(:,1);
    if size(info, 2) > 1 % if color col
      Classes.colors = info(:,2);
    end
    if size(info, 2) > 2 % if marker col
      Classes.markers = info(:,3);
    end
  catch
    Classes = struct();
    Classes.names = unique(classResults);
    Classes.colors = distinguishable_colors(length(Classes));
  end
  
  %% Transform Data
  data = struct();

  % Linear Data
  data.Linear.data = [simIDs classResults];
  
  variedNames = VariedData(1).varied(:);
  data.Linear.dimNames = [{'simID', 'class'}, variedNames(:)'];
  
  data.Linear.dimTypes = {'index', 'categorical'};
  data.Linear.dimTypes(end+1:end+length(variedNames)) = deal({'ordinal'});
  
  % Table Data
  simIDcharCell = cellfun(@num2str, simIDs,'uni',0);
  data.Table = table(categorical(classResults),'VariableNames',{'class'}, 'RowNames',simIDcharCell);
  data.Table.Properties.UserData.nonAxisDims = 1;
  
  for fld = variedNames(:)'
    if isnumeric(VariedData(1).(fld{1}))
      tempData = [VariedData.(fld{1})];
      data.Table.(fld{1}) = tempData(:);
      
      tempData = {VariedData.(fld{1})};
    else
      tempData = {VariedData.(fld{1})};
      data.Table.(fld{1}) = tempData(:);
    end
    data.Linear.data(:, end+1) = tempData(:);
  end

  % MultiDim Data
  ordinalDimInd = strcmp(data.Linear.dimTypes, 'ordinal');
  ordinalVars = data.Linear.dimNames(ordinalDimInd);
  nonOrdinalVars = data.Linear.dimNames(~ordinalDimInd);
  for iVar = 1:length(ordinalVars)
    thisVar = ordinalVars{iVar};
    varData{iVar} = data.Table.(thisVar);
    uniqueVals{iVar} = sort(unique(varData{iVar}));
    nVals(iVar) = length(uniqueVals{iVar});
  end
  
  matSize = cellfun(@length, uniqueVals);
  lin2mdMat = nan(matSize);
  
  for iData = 1:length(simIDs)
    mdInd = cell(1,length(ordinalVars));
    for iVar = 1:length(ordinalVars)
%       thisVar = ordinalVars{iVar};
      thisVarData = varData{iVar};
      mdInd{iVar} = find(uniqueVals{iVar} == thisVarData(iData));
    end
    lin2mdMat(mdInd{:})=iData;
  end
  lin2mdMat(isnan(lin2mdMat)) = length(simIDs)+1; %assign nan's to a new fake index
  
  for iVar = 1:length(nonOrdinalVars)
    thisNonOrdVar = nonOrdinalVars{iVar};
    thisLinearData = data.Linear.data(:, strcmp(data.Linear.dimNames,thisNonOrdVar));
    thisLinearData{end+1} = nan; %make fake index at end
    data.MultiDim.data{iVar} = thisLinearData(lin2mdMat);
    
    % Fix NaNs
    if strcmp(data.Linear.dimTypes{iVar}, 'index')
      % Assign -1 to index for missing values
      [data.MultiDim.data{iVar}{cellfun(@isnan, data.MultiDim.data{iVar})}] = deal(-1);
    elseif strcmp(data.Linear.dimTypes{iVar}, 'categorical')
      % Assign empty string for missing values
      sliceInd = cell2mat(cellfun(@isnanStr, data.MultiDim.data{iVar},'uni',0));
      [data.MultiDim.data{iVar}{sliceInd}] = deal('');
    end
  end
  
  
  data.MultiDim.dataNames = {'simID', 'class'};
  data.MultiDim.dataTypes = {'index', 'categorical'};
  data.MultiDim.dimNames = ordinalVars;
  data.MultiDim.dimVals = uniqueVals; % sorted
  data.MultiDim.nDimVals = nVals;
  
  % Labels
  Classes.linearDimNum = 2;
  Classes.multiDimNum = 2;
  Classes.tableName = 'class';
  data.Label = Classes;
  
  % Save
  save(filePath, 'data')
  vfprintf('\tSaved data.\n')
else % data file exists
  vfprintf('Loading data from file...\n')
  data = load(filePath);
  data = data.data;
  vfprintf('\tLoaded data from file.\n')
end


%% Sub Fns
  function output = isnanStr(x)
    if ischar(x) || isstring(x)
      output = false;
    else
      output = logical(isnan(x));
    end
  end

  function vfprintf(varargin)
    if options.verboseBool
      fprintf(varargin{:})
    end
  end
end
