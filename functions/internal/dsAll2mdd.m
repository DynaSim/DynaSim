function [xp,is_image] = dsAll2mdd(data,merge_covaried_axes,merge_sparse_axes,merge_everything,varargin)

if isempty(data)
  error('Input data is empty');
end

if nargin < 2
    merge_covaried_axes = true;
end

if nargin < 3
    merge_sparse_axes = true;
end

if nargin < 4
    merge_everything = false;
end

if ~isfield(data,'plot_files')      % Standard DynaSim data structure

    % Check inputs
    data=dsCheckData(data, varargin{:});
      % note: calling dsCheckData() at beginning enables analysis/plotting functions to
      % accept data matrix [time x cells] in addition to DynaSim data structure.

    % Convert input data to MDD
    xp = ds2mdd(data,merge_covaried_axes,merge_sparse_axes,merge_everything);
    is_image = 0;
else                            % Structure of links to plots

    % Convert input data to MDD
    data_img=data;
    xp = dsImg2mdd(data_img,merge_covaried_axes,merge_sparse_axes,merge_everything);
    
    is_image = 1;

end
