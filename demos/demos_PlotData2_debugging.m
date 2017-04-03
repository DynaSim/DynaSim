
%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % Advanced testing for debugging PlotData2 % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% Advanced testing for debugging - set up data

% Make sure data is loaded and path is ready
if ~exist('SimulateModel','file'); error('Use demos_PlotData2.m script to set up simulation path.'); end
if ~exist('data','var'); error('Use demos_PlotData2.m script to load data'); end

xp = DynaSim2xPlt(data);
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

data_single = xPlt2DynaSim(xp_single);
data_single_squeeze = xPlt2DynaSim(squeeze(xp_single));
data_single_pops = xPlt2DynaSim(xp_single_pops); data_single_pops=rmfield(data_single_pops,{'varied','I_E_tauD','E_Iapp'});  % This removes all the "vary" fields
data_single_vars = xPlt2DynaSim(xp_single_vars); data_single_vars=rmfield(data_single_vars,{'varied','I_E_tauD','E_Iapp'});  % This removes all the "vary" fields
data_row = xPlt2DynaSim(xp_row);
data_col = xPlt2DynaSim(xp_col);
data_col_pops = xPlt2DynaSim(xp_col_pops);
data_col_vars = xPlt2DynaSim(xp_col_vars);
data_col_varspops = xPlt2DynaSim(xp_col_varspops);
data_mat_pops = xPlt2DynaSim(xp_mat_pops);
data_mat_vars = xPlt2DynaSim(xp_mat_vars);
data_all = xPlt2DynaSim(xp);

vis = 'off';

%% Plot waveforms, comparing PlotData and PlotData2
close all; d = data_single; PlotData2(d,'visible',vis); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'visible',vis); %PlotData(d,'visible',vis);
close all; d = data_single_squeeze; PlotData2(d,'visible',vis); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'visible',vis); %PlotData(d,'visible',vis);
close all; d = data_single_pops; PlotData2(d,'visible',vis); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'visible',vis); PlotData(d,'visible',vis);
close all; d = data_single_vars; PlotData2(d,'visible',vis); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'visible',vis); %PlotData(d,'visible',vis);
close all; d = data_row; PlotData2(d,'visible',vis); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'visible',vis); PlotData(d,'visible',vis);
close all; d = data_col; PlotData2(d,'visible',vis); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'visible',vis); PlotData(d,'visible',vis);
close all; d = data_col_pops; PlotData2(d,'visible',vis); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'visible',vis); PlotData(d,'visible',vis);
close all; d = data_col_vars; PlotData2(d,'visible',vis); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'visible',vis); PlotData(d,'visible',vis);
close all; d = data_col_varspops; PlotData2(d,'visible',vis); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'visible',vis); PlotData(d,'visible',vis);
close all; d = data_mat_pops; PlotData2(d,'do_mean',1); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'do_mean',1); PlotData(d,'visible',vis);
close all; d = data_mat_vars; PlotData2(d,'visible',vis); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'visible',vis); PlotData(d,'visible',vis);
close all; d = data_all; PlotData2(d,'do_mean',1); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'do_mean',1); PlotData(d,'visible',vis);


%% Plot rastergrams, comparing PlotData and PlotData2
close all; d = data_single_pops; PlotData2(d,'plot_type','rastergram','visible',vis); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'plot_type','rastergram','visible',vis); PlotData(d,'plot_type','rastergram','visible',vis);
close all; d = data_single_vars; PlotData2(d,'plot_type','rastergram','visible',vis); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'plot_type','rastergram','visible',vis); PlotData(d,'plot_type','rastergram','visible',vis);
close all; d = data_row; PlotData2(d,'plot_type','rastergram','visible',vis); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'plot_type','rastergram','visible',vis); PlotData(d,'plot_type','rastergram','visible',vis);
close all; d = data_col; PlotData2(d,'plot_type','rastergram','visible',vis); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'plot_type','rastergram','visible',vis); PlotData(d,'plot_type','rastergram','visible',vis);
close all; d = data_col_pops; PlotData2(d,'plot_type','rastergram','visible',vis); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'plot_type','rastergram','visible',vis); PlotData(d,'plot_type','rastergram','visible',vis);
close all; d = data_col_vars; PlotData2(d,'plot_type','rastergram','visible',vis); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'plot_type','rastergram','visible',vis); PlotData(d,'plot_type','rastergram','visible',vis);
close all; d = data_col_varspops; PlotData2(d,'plot_type','rastergram','visible',vis); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'plot_type','rastergram','visible',vis); PlotData(d,'plot_type','rastergram','visible',vis);
% close all; d = data_mat_pops; PlotData2(d,'plot_type','rastergram','do_mean',1); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'plot_type','rastergram','do_mean',1); PlotData(d,'plot_type','rastergram','visible',vis);
close all; d = data_mat_vars; PlotData2(d,'plot_type','rastergram','visible',vis); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'plot_type','rastergram','visible',vis); PlotData(d,'plot_type','rastergram','visible',vis);
% close all; d = data_all; PlotData2(d,'plot_type','rastergram','do_mean',1); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'plot_type','rastergram','do_mean',1); PlotData(d,'plot_type','rastergram','visible',vis);


%% Plot FR2, comparing PlotData and PlotData2
close all; d = data_single_pops; PlotData2(d,'plot_type','heatmap_sortedFR','visible',vis); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'plot_type','heatmap_sortedFR','visible',vis); PlotFR2(d,'plot_type','heatmap_sorted','visible',vis);
close all; d = data_single_vars; PlotData2(d,'plot_type','heatmap_sortedFR','visible',vis); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'plot_type','heatmap_sortedFR','visible',vis); PlotFR2(d,'plot_type','heatmap_sorted','visible',vis);
close all; d = data_row; PlotData2(d,'plot_type','heatmap_sortedFR','visible',vis); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'plot_type','heatmap_sortedFR','visible',vis); PlotFR2(d,'plot_type','heatmap_sorted','visible',vis);
close all; d = data_col; PlotData2(d,'plot_type','heatmap_sortedFR','visible',vis); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'plot_type','heatmap_sortedFR','visible',vis); PlotFR2(d,'plot_type','heatmap_sorted','visible',vis);
close all; d = data_col_pops; PlotData2(d,'plot_type','heatmap_sortedFR','visible',vis); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'plot_type','heatmap_sortedFR','visible',vis); PlotFR2(d,'plot_type','heatmap_sorted','visible',vis);
close all; d = data_col_vars; PlotData2(d,'plot_type','heatmap_sortedFR','visible',vis); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'plot_type','heatmap_sortedFR','visible',vis); PlotFR2(d,'plot_type','heatmap_sorted','visible',vis);
close all; d = data_col_varspops; PlotData2(d,'plot_type','heatmap_sortedFR','visible',vis); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'plot_type','heatmap_sortedFR','visible',vis); PlotFR2(d,'plot_type','heatmap_sorted','visible',vis);
% close all; d = data_mat_pops; PlotData2(d,'plot_type','heatmap_sortedFR','do_mean',1); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'plot_type','heatmap_sortedFR','do_mean',1); PlotFR2(d,'plot_type','heatmap_sorted','visible',vis);
close all; d = data_mat_vars; PlotData2(d,'plot_type','heatmap_sortedFR','visible',vis); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'plot_type','heatmap_sortedFR','visible',vis); PlotFR2(d,'plot_type','heatmap_sorted','visible',vis);
% close all; d = data_all; PlotData2(d,'plot_type','heatmap_sortedFR','do_mean',1); d = xPlt2DynaSim(DynaSim2xPlt(d)); PlotData2(d,'plot_type','heatmap_sortedFR','do_mean',1); PlotFR2(d,'plot_type','heatmap_sorted','visible',vis);



keyboard

%% Advanced testing for debugging - set up data_img

xp_img = DynaSimImg2xPlt(data_img);
