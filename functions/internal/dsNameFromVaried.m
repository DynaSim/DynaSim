function new_filename = dsNameFromVaried(data, prefix, old_filename)
%NAMEFROMVARIED - makes a filename based on the parameters in data.varied.
%
% Usage:
%   new_filename = dsNameFromVaried(data, file_type, old_filename)
%
% Inputs:
%   - data: DynaSim data structure (also accepted: data file name)
%   - prefix: file prefix for type of file
%   - result_file: previous result_file, only uses this for the path.
%
% Outputs:
%   - new_filename: where to save result, based on file_type and data.varied

pathstr = fileparts2(old_filename);

fileName = prefix;

%check for simID# from batch sims
token = regexp(old_filename, '(sim\d+)', 'tokens');
if ~isempty(token)
  fileName = [fileName '_' token{1}{1}];
end

  for param = data.varied(:)'
    fileName = [fileName '__' param{1} '_' sprintf('%g',data.(param{1}))];
  end

new_filename = fullfile(pathstr, fileName);

end
