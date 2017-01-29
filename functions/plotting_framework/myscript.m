
% This is a script that I run through every time I want to test the code.
% Assumes I've just run one of the scripts in demos.m, namely the section
% titled "%% Save data from a set of simulations"
% Then, I load demo_sPING_3:
% data=ImportData('demo_sPING_3');


%% Load the data and import into xPlt class
% ...Assumes we have some DynaSim data already loaded...
% data=ImportData('demo_sPING_3');
% cd ../../functions/plotting_framework

% Load the data linearly
[data_linear,ax,ax_names,time] = DynaSimExtract (data);



% Import into an xPlt class
xp = xPlt;
xp = xp.importLinearData(data_linear,ax{:});
xp = xp.importAxisNames(ax_names);


%% Try selecting a subset
xp2 = xp.subset(2,2,[],7:8);

xp3 = xp;
xp3.data = xp.data(:,:,:,5:8);      % Brute force selection of the xp.data.    
                                    % Should produce an error when run
                                    % xp3.squeeze since axes dimensions
                                    % mismatch.

%% Try selecting another subset for actual plotting
clear xp2 xp3
xp4 = squeeze(xp.subset(3,[],[],8));
recursivePlot(xp4,{@xp_handles_subplotrows,@xp_images_subplotcols,@xp_matrix_basicplot},{1,1,1});



