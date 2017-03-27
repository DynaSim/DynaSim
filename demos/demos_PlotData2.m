%% % % % % % % % PlotData2 tutorial % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%% Set up paths 
% Get ready...

% Add DynaSim to path if it's not already there
if exist('setup_DynaSim_path','file')
    setup_DynaSim_path;
else
    error('Add the DynaSim folder to the MATLAB path - e.g. run addpath(genpath(DynaSimPath))');
end

% Set where to save outputs
output_directory = getDsVar('demos_path');
study_dir = fullfile(output_directory,'demo_sPING_100cells_3x3');

%% Load data
% Load data in traditional DynaSim format
data=ImportData(study_dir);


%% Do some basic plots with PlotData2

% Default settings
PlotData2(data);

% Take mean
PlotData2(data,'do_mean',1);

% Plot data, just E cells, with zoom turned on
PlotData2(data,'do_zoom',1,'population','E');

% As above, but focus on only a subset of the varied data. Also, reduce the
% total number of traces displayed.
PlotData2(data,'population','E','variable','v','E_Iapp',1:3,'varied2',[1,3],'do_zoom',1,'max_num_overlaid',3);

%% Recursive plots with PlotData2

% The number of subplots packed into a single figure can be substantially
% increased by recurisve subplotting. This is controlled with the
% "num_embedded_subplots" flag, which can range from 1-4.

% Plot membrane voltage for E and I cells across the parameter sweep
% % (1 variable, 2 pops, varied1, varied2: Ndims = 3)
PlotData2(data,'num_embedded_subplots',2,'do_zoom',1);  % Default
PlotData2(data,'num_embedded_subplots',4,'do_zoom',1);  % Increase embedding

% Note that in the prev example, although we set embedded subplots to 4, only
% the first 3 are used. If excess subplots are requested and unused, they
% will be thrown out.

% Setting num_embedded_subplots to provides the same data but in a different
% arrangement.
PlotData2(data,'num_embedded_subplots',3,'do_zoom',1);

% Depending on the data being plotted, different some arrangements can be
% more useful than others. For example, this plot embeds varied as a
% sub-subplot.
% % (2 variables, 2 pops, varied2, (varied1 is fixed): Ndims = 3)
PlotData2(data,'num_embedded_subplots',3,'variable','iNa*','varied1',2)

% Lastly, 4 subplots can be nested together to view 4 dimensions
% simultaneously. This can be increased to 5 with the overlay function
% described below, and 6 if you count multiple figures. While these figures
% may get crowded, the supersize_me flag, also described below, aims to 
% address this by producing single figures occupying very large canvases,
% which can be zoomed to a high level of detail.
% % (2 vars, 2 pops, varied1, varied2: Ndims = 4)
PlotData2(data,'max_num_overlaid',3,'num_embedded_subplots',4,'population','all','variable','v|iNa_h','varied1',2:3,'lock_axes',false);


%% Modifying overlaid traces

% This provides an alternative method for stacking more data into a single
% plot.

% If there is only 1 cell from each population, PlotData will overlay
% cells from different populations by default.
% This can be achieved by averaging across cells...
PlotData2(data,'do_mean',1,'varied1',1:3,'varied2',2:3)

% ... or by looking at only the 1st cell in each population
PlotData2(data,'max_num_overlaid',1,'varied1',1:3,'varied2',2:3)

% Turning this overlay off splits up the figures
PlotData2(data,'max_num_overlaid',1,'varied1',1:3,'varied2',2:3,'force_overlay','none')

% Overlays can also be doubly stacked
PlotData2(data,'max_num_overlaid',10,'force_overlay','populations');

% However, these can be messy, so the groups can be shifted up or down
PlotData2(data,'max_num_overlaid',10,'force_overlay','populations','do_overlay_shift',true);
    % Variables are stacked from top to bottom, so E cells are on top and I
    % cells underneath

% You can overlay any variable. For example, here is injected current
PlotData2(data,'do_mean',0,'varied1',1:3,'varied2',2:3,'force_overlay','varied1','do_overlay_shift',1,'overlay_shift_val',100)


% Variables with very different units can be compared side-by-side by
% taking the z-score first. In this case we compare E cell membrane voltage
% to its inhibitory synaptic input.
PlotData2(data,'population','E','variable','v|I_iGABAa_s','force_overlay','variables','do_overlay_shift',true,'overlay_shift_val',3,'do_zscore',true,'do_zoom',1);

% % Double stack overlays
% PlotData2(data,'population','E','variable','iNa_m|I_iGABAa_s','force_overlay','variables','do_overlay_shift',true);
% PlotData2(data,'population','I','variable','v');


%% Import pre-saved images

% Import pre-saved images and tile them
PlotData2(study_dir)

% As above, but supersize them
PlotData2(study_dir,'supersize_me',1)


%% Supersize me
load sample_data_dynasim_large.mat

PlotData2(data);
