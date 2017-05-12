function [xp,is_image] = all2MDD(data,varargin)

if ~isfield(data,'plot_files')      % Standard DynaSim data structure

    % Check inputs
    data=ds.checkData(data, varargin{:});
      % note: calling ds.checkData() at beginning enables analysis/plotting functions to
      % accept data matrix [time x cells] in addition to DynaSim data structure.


    % Convert input data to MDD
    xp = ds.ds2MDD(data);
    is_image = 0;
else                            % Structure of links to plots

    % Convert input data to MDD
    data_img=data;
    xp = ds.img2MDD(data_img);
    
    is_image = 1;

end
