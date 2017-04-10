
%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % Advanced testing for debugging ds.plotData2 % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% Advanced testing for debugging - set up data

% Make sure data is loaded and path is ready
if ~exist('ds.simulateModel','file'); error('Use demos_ds.plotData2.m script to set up simulation path.'); end
if ~exist('data','var'); error('Use demos_ds.plotData2.m script to load data'); end

xp = ds.ds2xPlt(data);
xp_single = xp(1,1,1,'v');
xp_single_pops = xp(1,1,:,'v');
xp_single_vars = xp(1,1,1,:);
xp_row = xp(:,1,1,'v');
xp_col = xp(1,:,1,'v');
xp_col_pops = xp(2,:,:,'v');
xp_col_vars = xp(1,:,1,:);
xp_col_varspops = xp(1,:,:,:);
xp_mat_pops = xp(:,:,:,'v');
xp_mat_vars = xp(:,:,1,:);

data_single = ds.xPlt.xPlt2ds(xp_single);
data_single_squeeze = ds.xPlt.xPlt2ds(squeeze(xp_single));
data_single_pops = ds.xPlt.xPlt2ds(xp_single_pops); data_single_pops=rmfield(data_single_pops,{'varied','I_E_tauD','E_Iapp'});  % This removes all the "vary" fields
data_single_vars = ds.xPlt.xPlt2ds(xp_single_vars); data_single_vars=rmfield(data_single_vars,{'varied','I_E_tauD','E_Iapp'});  % This removes all the "vary" fields
data_row = ds.xPlt.xPlt2ds(xp_row);
data_col = ds.xPlt.xPlt2ds(xp_col);
data_col_pops = ds.xPlt.xPlt2ds(xp_col_pops);
data_col_vars = ds.xPlt.xPlt2ds(xp_col_vars);
data_col_varspops = ds.xPlt.xPlt2ds(xp_col_varspops);
data_mat_pops = ds.xPlt.xPlt2ds(xp_mat_pops);
data_mat_vars = ds.xPlt.xPlt2ds(xp_mat_vars);
data_all = ds.xPlt.xPlt2ds(xp);

vis = 'off';

%% Plot waveforms, comparing ds.plotData and ds.plotData2
close all; d = data_single; ds.plotData2(d,'visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'visible',vis); %ds.plotData(d,'visible',vis);
close all; d = data_single_squeeze; ds.plotData2(d,'visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'visible',vis); %ds.plotData(d,'visible',vis);
close all; d = data_single_pops; ds.plotData2(d,'visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'visible',vis); ds.plotData(d,'visible',vis);
close all; d = data_single_vars; ds.plotData2(d,'visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'visible',vis); %ds.plotData(d,'visible',vis);
close all; d = data_row; ds.plotData2(d,'visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'visible',vis); ds.plotData(d,'visible',vis);
close all; d = data_col; ds.plotData2(d,'visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'visible',vis); ds.plotData(d,'visible',vis);
close all; d = data_col_pops; ds.plotData2(d,'visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'visible',vis); ds.plotData(d,'visible',vis);
close all; d = data_col_vars; ds.plotData2(d,'visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'visible',vis); ds.plotData(d,'visible',vis);
close all; d = data_col_varspops; ds.plotData2(d,'visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'visible',vis); ds.plotData(d,'visible',vis);
close all; d = data_mat_pops; ds.plotData2(d,'do_mean',1,'visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'do_mean',1,'visible',vis); ds.plotData(d,'visible',vis);
close all; d = data_mat_vars; ds.plotData2(d,'visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'visible',vis); ds.plotData(d,'visible',vis);
close all; d = data_all; ds.plotData2(d,'do_mean',1,'visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'do_mean',1,'visible',vis); ds.plotData(d,'visible',vis);


%% Plot rastergrams, comparing ds.plotData and ds.plotData2
close all; d = data_single_pops; ds.plotData2(d,'plot_type','rastergram','visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'plot_type','rastergram','visible',vis); ds.plotData(d,'plot_type','rastergram','visible',vis);
close all; d = data_single_vars; ds.plotData2(d,'plot_type','rastergram','visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'plot_type','rastergram','visible',vis); ds.plotData(d,'plot_type','rastergram','visible',vis);
close all; d = data_row; ds.plotData2(d,'plot_type','rastergram','visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'plot_type','rastergram','visible',vis); ds.plotData(d,'plot_type','rastergram','visible',vis);
close all; d = data_col; ds.plotData2(d,'plot_type','rastergram','visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'plot_type','rastergram','visible',vis); ds.plotData(d,'plot_type','rastergram','visible',vis);
close all; d = data_col_pops; ds.plotData2(d,'plot_type','rastergram','visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'plot_type','rastergram','visible',vis); ds.plotData(d,'plot_type','rastergram','visible',vis);
close all; d = data_col_vars; ds.plotData2(d,'plot_type','rastergram','visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'plot_type','rastergram','visible',vis); ds.plotData(d,'plot_type','rastergram','visible',vis);
close all; d = data_col_varspops; ds.plotData2(d,'plot_type','rastergram','visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'plot_type','rastergram','visible',vis); ds.plotData(d,'plot_type','rastergram','visible',vis);
% close all; d = data_mat_pops; ds.plotData2(d,'plot_type','rastergram','do_mean',1,'visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'plot_type','rastergram','do_mean',1,'visible',vis); ds.plotData(d,'plot_type','rastergram','visible',vis);
close all; d = data_mat_vars; ds.plotData2(d,'plot_type','rastergram','visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'plot_type','rastergram','visible',vis); ds.plotData(d,'plot_type','rastergram','visible',vis);
% close all; d = data_all; ds.plotData2(d,'plot_type','rastergram','do_mean',1,'visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'plot_type','rastergram','do_mean',1,'visible',vis); ds.plotData(d,'plot_type','rastergram','visible',vis);


%% Plot FR2, comparing ds.plotData and ds.plotData2
close all; d = data_single_pops; ds.plotData2(d,'plot_type','heatmap_sortedFR','visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'plot_type','heatmap_sortedFR','visible',vis); ds.plotFR2(d,'plot_type','heatmap_sorted','visible',vis);
close all; d = data_single_vars; ds.plotData2(d,'plot_type','heatmap_sortedFR','visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'plot_type','heatmap_sortedFR','visible',vis); ds.plotFR2(d,'plot_type','heatmap_sorted','visible',vis);
close all; d = data_row; ds.plotData2(d,'plot_type','heatmap_sortedFR','visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'plot_type','heatmap_sortedFR','visible',vis); ds.plotFR2(d,'plot_type','heatmap_sorted','visible',vis);
close all; d = data_col; ds.plotData2(d,'plot_type','heatmap_sortedFR','visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'plot_type','heatmap_sortedFR','visible',vis); ds.plotFR2(d,'plot_type','heatmap_sorted','visible',vis);
close all; d = data_col_pops; ds.plotData2(d,'plot_type','heatmap_sortedFR','visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'plot_type','heatmap_sortedFR','visible',vis); ds.plotFR2(d,'plot_type','heatmap_sorted','visible',vis);
close all; d = data_col_vars; ds.plotData2(d,'plot_type','heatmap_sortedFR','visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'plot_type','heatmap_sortedFR','visible',vis); ds.plotFR2(d,'plot_type','heatmap_sorted','visible',vis);
close all; d = data_col_varspops; ds.plotData2(d,'plot_type','heatmap_sortedFR','visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'plot_type','heatmap_sortedFR','visible',vis); ds.plotFR2(d,'plot_type','heatmap_sorted','visible',vis);
% close all; d = data_mat_pops; ds.plotData2(d,'plot_type','heatmap_sortedFR','do_mean',1,'visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'plot_type','heatmap_sortedFR','do_mean',1,'visible',vis); ds.plotFR2(d,'plot_type','heatmap_sorted','visible',vis);
close all; d = data_mat_vars; ds.plotData2(d,'plot_type','heatmap_sortedFR','visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'plot_type','heatmap_sortedFR','visible',vis); ds.plotFR2(d,'plot_type','heatmap_sorted','visible',vis);
% close all; d = data_all; ds.plotData2(d,'plot_type','heatmap_sortedFR','do_mean',1,'visible',vis); d = ds.xPlt.xPlt2ds(ds.ds2xPlt(d)); ds.plotData2(d,'plot_type','heatmap_sortedFR','do_mean',1,'visible',vis); ds.plotFR2(d,'plot_type','heatmap_sorted','visible',vis);



keyboard

%% Advanced testing for debugging - set up data_img

xp_img = ds.img2xPlt(data_img);
