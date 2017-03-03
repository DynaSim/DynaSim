function result = SaveVaried(data, varargin)
%% SaveVaried
% Author: Erik Roberts
%
% Purpose: Saves the data fields listed in varied.

if numel(data)>1
  % use AnalyzeStudy to recursively call ClassifyTCRE on each data set
  result=AnalyzeStudy(data,@SaveVaried,varargin{:});
  return;
end

fields = data.varied;

result = struct();
for field = fields(:)'
  field = field{1};
  result.(field) = data.(field);
end

end