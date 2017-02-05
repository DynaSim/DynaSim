
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


%% Try selecting a bunch of different subsets, squeezing, and testing for errors
% This should catch most bugs that crop up.
xp2 = xp.subset(2,2,[],7:8);    % 1x1x2x2
    xp2.checkDims;
    xp2.getaxisinfo
    xp2 = xp2.squeeze;
    xp2.getaxisinfo
    xp2.checkDims
    %%
clc
xp2 = xp.subset(2,2,[],8);      % 1x1x2x1
    xp2.checkDims;
    xp2.getaxisinfo
    xp2 = xp2.squeeze;
    xp2.getaxisinfo
    xp2.checkDims
    %%
clc
xp2 = xp.subset(2,2,2,[]);      % 1x1x1x8
    xp2.checkDims;
    xp2.getaxisinfo
    xp2 = xp2.squeeze;
    xp2.getaxisinfo
    xp2.checkDims
    %%
clc
xp2 = xp.subset(1,[],2,1);      % 1x2x1x1
    xp2.checkDims;
    xp2.getaxisinfo
    xp2 = xp2.squeeze;
    xp2.getaxisinfo
    xp2.checkDims
    %%
clc
xp2 = xp.subset([],1,2,1);      % 2x1x1x1
    xp2.checkDims;
    xp2.getaxisinfo
    xp2 = xp2.squeeze;
    xp2.getaxisinfo
    xp2.checkDims
 
    %%
clc
xp3 = xp;
xp3.data = xp.data(:,:,:,5:8);      % Brute force selection of the xp.data.    
                                    % Should produce an error when run
                                    % xp3.squeeze since axes dimensions
                                    % mismatch.
xp3.getaxisinfo;
% xp3.checkDims;      % Should produce an error
xp3 = xp3.fixAxes;
xp3.checkDims;
xp3.getaxisinfo;
                                    
xp6 = xp.subset(1,[],2,1);

%% Try selecting another subset for actual plotting
clear xp2 xp3
xp4 = (xp.subset([],[],[],8));
%%
% recursivePlot(xp4,{@xp_subplot,@xp_subplot,@xp_matrix_basicplot},{1:2,3},{{[],1},{1,1},{}});
% recursivePlot(xp4,{@xp_subplot_grid3D,@xp_subplot,@xp_matrix_basicplot},{1:2,3},{{},{0,1},{}});
% recursivePlot(xp4,{@xp_subplot_grid3D,@xp_matrix_basicplot},{[3,1,2]},{{},{}});
recursivePlot(xp4,{@xp_subplot_grid3D,@xp_subplot,@xp_matrix_basicplot},{[1,2,4],3},{{},{0,1},{}});



%% Test subset with regular expressions
xp5 = xp.subset([],[],[1],'iNa*');

%% Test packDims
clear xp2 xp3 xp4 xp5
xp2 = xp.subset(2,2,[],[1,3,5:8]);
xp2 = xp2.squeeze;
%%
xp2 = xp2.packDim(2,3);

