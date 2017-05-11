function [xp,is_image] = all2MDict(data,varargin)

if ~isfield(data,'plot_files')      % Standard DynaSim data structure

    % Check inputs
    data=ds.checkData(data, varargin{:});
      % note: calling ds.checkData() at beginning enables analysis/plotting functions to
      % accept data matrix [time x cells] in addition to DynaSim data structure.


    % Convert input data to MDict
    xp = ds.ds2MDict(data);
    is_image = 0;
else                            % Structure of links to plots

    % Convert input data to MDict
    data_img=data;
    xp = ds.img2MDict(data_img);
    
    is_image = 1;

end
