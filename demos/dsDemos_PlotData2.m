%% % % % % % % % dsPlot2 tutorial % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%% Set up paths 
% Get ready...

% Add DynaSim to path if it's not already there
if exist('setupDynaSimPath','file')
    setupDynaSimPath;
else
    error('Add the DynaSim folder to the MATLAB path - e.g. run addpath(genpath(DynaSimPath))');
end

% Set where to save outputs
output_directory = ds.getConfig('demos_path');
study_dir = fullfile(output_directory,'demo_sPING_100cells_3x3');

%% Import the data

% Make sure sample data exists; if not copy it into place
if ~exist(study_dir,'dir')
    ds.unzipDemoData(study_dir);
end

% Load data in traditional DynaSim format
data=ds.importData(study_dir);

% % Import saved plots from huge sweep of images
% data_3D_plots = ds.importPlots(fullfile(output_directory,'demo_sPING_3b_3D'));


%% Do some basic plots with dsPlot2

close all

% Default settings
dsPlot2(data);

% Take mean
dsPlot2(data,'do_mean',1);

% Plot data, just E cells, with zoom turned on
dsPlot2(data,'do_zoom',1,'population','E');

% As above, but focus on only a subset of the varied data. Also, reduce the
% total number of traces displayed.
dsPlot2(data,'population','E','variable','v','E_Iapp',1:3,'varied2',[1,3],'do_zoom',1,'max_num_overlaid',10);


% All basic plot_types from dsPlot should work
% Options are availalbe from dsPlot:
%     {'waveform','imagesc','rastergram','raster','power'}
% And also from ds.plotFR2
%     {'heatmapFR','heatmap_sortedFR','meanFR','meanFRdens'}
dsPlot2(data,'population','all','plot_type','rastergram')              % Rastergram
dsPlot2(data,'population','I','plot_type','heatmap_sortedFR')          % Firing rate (FR) heatmap

%% Recursive plots with dsPlot2

close all

% The number of subplots packed into a single figure can be substantially
% increased by recurisve subplotting. This is controlled with the
% "num_embedded_subplots" flag, which can range from 1-4.

% Plot membrane voltage for E and I cells across the parameter sweep
% % (1 variable, 2 pops, varied1, varied2: Ndims = 3)
dsPlot2(data,'num_embedded_subplots',2,'do_zoom',1,'max_num_overlaid',10);  % Default
dsPlot2(data,'num_embedded_subplots',4,'do_zoom',1,'max_num_overlaid',10);  % Nested embedding

% Note that in the prev example, although we set embedded subplots to 4, only
% the first 3 are used. If excess subplots are requested and unused, they
% will be thrown out.
% 
% % Setting num_embedded_subplots to provides the same data but in a different
% % arrangement.
% dsPlot2(data,'num_embedded_subplots',3,'do_zoom',1);
% 
% % Depending on the data being plotted, different some arrangements can be
% % more useful than others. For example, this plot embeds varied as a
% % sub-subplot.
% % % (2 variables, 2 pops, varied2, (varied1 is fixed): Ndims = 3)
% dsPlot2(data,'num_embedded_subplots',3,'variable','iNa*','varied1',2)
% 
% Lastly, 4 subplots can be nested together to view 4 dimensions
% simultaneously. This can be increased to 5 with the overlay function
% described below, and 6 if you count multiple figures. While these figures
% may get crowded, the supersize_me flag, also described below, aims to 
% address this by producing single figures occupying very large canvases,
% which can be zoomed to a high level of detail.
% % (2 vars, 2 pops, varied1, varied2: Ndims = 4)
dsPlot2(data,'max_num_overlaid',3,'num_embedded_subplots',4,'population','all','variable','v|iNa_h','varied1',2:3,'lock_axes',false);


%% Modifying overlaid traces

close all

% This provides an alternative method for stacking more data into a single
% plot.

% If there is only 1 cell from each population, dsPlot will overlay
% cells from different populations by default.
% This can be achieved by averaging across cells...
dsPlot2(data,'do_mean',1,'varied1',1:3,'varied2',2:3)

% ... or by looking at only the 1st cell in each population
dsPlot2(data,'max_num_overlaid',1,'varied1',1:3,'varied2',2:3)

% Turning this overlay off splits up the figures
dsPlot2(data,'max_num_overlaid',1,'varied1',1:3,'varied2',2:3,'force_overlay','none')

close all;

% Overlays can also be doubly stacked. Note that both E and I cells are
% shown on a single plot.
dsPlot2(data,'max_num_overlaid',5,'force_overlay','populations');

% However, these can be messy, so the groups can be shifted up or down
dsPlot2(data,'max_num_overlaid',5,'force_overlay','populations','do_overlay_shift',true);
    % Variables are stacked from top to bottom, so E cells are on top and I
    % cells underneath

% You can overlay any variable. For example, here is injected current.
% This is accomplished by setting the axis "varied1", which corresponds to
% the injected current, to be the variable that is rastered across the
% overlays.
dsPlot2(data,'do_mean',0,'varied1',1:3,'varied2',2:3,...
    'force_overlay','varied1','do_overlay_shift',1,'overlay_shift_val',100);


% Variables with very different units can be compared side-by-side by
% taking the z-score first. In this case we compare E cell membrane voltage
% to its inhibitory synaptic input.
dsPlot2(data,'population','E','variable','v|I_iGABAa_s','force_overlay','variables','do_overlay_shift',true,'overlay_shift_val',3,'do_zscore',true,'do_zoom',1);

% % Double stack overlays
% dsPlot2(data,'population','E','variable','iNa_m|I_iGABAa_s','force_overlay','variables','do_overlay_shift',true);
% dsPlot2(data,'population','I','variable','v');


%% Import pre-saved images
close all

% Import plot files
data_img = ds.importPlots(study_dir);
dsPlot2(data_img);

% Can also reference study_dir directly
dsPlot2(study_dir)

% As above, but supersize them
dsPlot2(study_dir,'supersize_me',1)


%% Supersize me
% load sample_data_dynasim_large.mat
% 
% dsPlot2(data);


%% Merging simulation output data

study_dir2 = fullfile(output_directory,'demo_sPING_100cells_5x1');

% Make sure sample data exists; if not copy it into place
if ~exist(study_dir2,'dir')
    ds.unzipDemoData(study_dir2);
end

% Load data in traditional DynaSim format
data_sim2=ds.importData(study_dir2);

% Test plot of new data
dsPlot2(data_sim2);

% Plot the merged new data with the old data
data_merged = ds.mergeData(data,data_sim2);
dsPlot2(data_merged,'do_mean',1)


%% Merge simulations' saved images

data_img_sim2=ds.importPlots(study_dir2);

% Test plot of saved images from sim2
dsPlot2(data_img_sim2);

% Merge
data_img_merged = ds.mergeData(data_img,data_img_sim2);
dsPlot2(data_img_merged);




