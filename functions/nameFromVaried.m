function new_result_file = nameFromVaried(data, plot_type, old_result_file)
%% new_result_file = nameFromVaried(data, plot_type, old_result_file)
% purpose: makes a filename based on the parameters in data.varied.
%
% inputs:
%   data: DynaSim data structure (also accepted: data file name)
%   plot_type: what to plot
%   result_file: previous result_file, only uses this for the path.
%
% outputs:
%   new_result_file: where to save result, based on plot_type and data.varied

pathstr = fileparts(old_result_file);

fileName = plot_type;
for param = data.varied(:)'
  fileName = [fileName '__' param{1} '_' num2str(data.(param{1}))];
end

new_result_file = fullfile(pathstr, fileName);

end