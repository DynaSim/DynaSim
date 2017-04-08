function [xp,is_image] = all2xPlt(data)

if ~isfield(data,'plot_files')      % Standard DynaSim data structure

    % Check inputs
    data=checkData(data);
      % note: calling checkData() at beginning enables analysis/plotting functions to
      % accept data matrix [time x cells] in addition to DynaSim data structure.


    % Convert input data to xPlt
    xp = DynaSim2xPlt(data);
    is_image = 0;
else                            % Structure of links to plots

    % Convert input data to xPlt
    data_img=data;
    xp = DynaSimImg2xPlt(data_img);
    
    is_image = 1;

end