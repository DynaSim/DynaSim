%% Set up paths 
% Get ready...

% Format
format compact

% Set path to your copy of the DynaSim toolbox
dynasim_path = fullfile(pwd, 'DynaSim')

% add DynaSim toolbox to Matlab path
addpath(genpath(dynasim_path)); % comment this out if already in path

% Set where to save outputs
output_directory = 'demos/outputs';

% move to root directory where outputs will be saved
output_fulldir = fullfile(dynasim_path, output_directory)
mkdir(output_fulldir);
cd(output_fulldir);

%% Run simulation - Sparse Pyramidal-Interneuron-Network-Gamma (sPING)

% define equations of cell model (same for E and I populations)
eqns={ 
  'dv/dt=Iapp+@current+noise*randn(1,N_pop)';
  'monitor iGABAa.functions, iAMPA.functions'
};
% Tip: monitor all functions of a mechanism using: monitor MECHANISM.functions

% create DynaSim specification structure
s=[];
s.populations(1).name='E';
s.populations(1).size=80;
s.populations(1).equations=eqns;
s.populations(1).mechanism_list={'iNa','iK'};
s.populations(1).parameters={'Iapp',5,'gNa',120,'gK',36,'noise',40};
s.populations(2).name='I';
s.populations(2).size=20;
s.populations(2).equations=eqns;
s.populations(2).mechanism_list={'iNa','iK'};
s.populations(2).parameters={'Iapp',0,'gNa',120,'gK',36,'noise',40};
s.connections(1).direction='I->E';
s.connections(1).mechanism_list={'iGABAa'};
s.connections(1).parameters={'tauD',10,'gSYN',.1,'netcon','ones(N_pre,N_post)'};
s.connections(2).direction='E->I';
s.connections(2).mechanism_list={'iAMPA'};
s.connections(2).parameters={'tauD',2,'gSYN',.1,'netcon',ones(80,20)};

% Vary two parameters (run a simulation for all combinations of values)
vary={
  'E'   ,'Iapp',[0 10 20];      % amplitude of tonic input to E-cells
  'I->E','tauD',[5 10 15]       % inhibition decay time constant from I to E
  };
SimulateModel(s,'save_data_flag',1,'study_dir','demo_sPING_3b',...
                'vary',vary,'verbose_flag',1, ...
                'save_results_flag',1,'plot_functions',@PlotData,'plot_options',{'format','png'} );



%% Load the data and import it into xPlt class

% Load data in traditional DynaSim format
data=ImportData('demo_sPING_3b');

% Extract the data in a linear table format
[data_table,column_titles,time] = Data2Table (data);

% Preview the contents of this table
%     Note: We cannot make this one big cell array since we want to allow
%     axis labels to be either strings or numerics.
previewTable(data_table,column_titles);

% Import the linear data into an xPlt object
xp = xPlt;
X = data_table{1};                          % X holds the data that will populate the multidimensional array. Must be numeric or cell array.
axislabels = data_table(2:end);             % Each entry in X has an associated set of axis labels, which will define its location in multidimensional space. **Must be numeric or cell array of chars only**
xp = xp.importLinearData(X,axislabels{:});
xp = xp.importAxisNames(column_titles(2:end));  % There should be 1 axis name for every axis, of type char.


% xPlt objects are essentially cell arrays (or matricies), but with the
% option to index using strings instead of just integers. 
% Thus, they are analogous to dictionaries in Python.
% (This core functionality is implemented by the multidimensional
% dictionaries (nDDict), which xPlt inherits, and to which xPlt adds
% plotting functionality.)
disp(xp);


% At its core, xPlt has 3 fields. xp.data stores the actual data (either a 
% matrix or a cell array). 
disp(xp.data);
size(xp.data)

% Next, xp.axis stores axis labels associated with each dimension in
% xp.data.
disp(xp.axis(1));

% Axis.values stores the actual axis labels. These can be numeric...
disp(xp.axis(1).values);

% ...or string type. As we shall see below, these axis labels can be
% referenced via index or regular expression.
disp(xp.axis(4).values);

% Axis.name field stores the name of the dimension. In this case, it is
% named after the parameter in the model that was varied.
disp(xp.axis(1).name)

% Axis.astruct is for internal use.
xp.axis(1).astruct

% xp.meta stores meta data for use by the user as they see fit.
% Here we will add some custom info to xp.metadata. This can be whatever
% you want. Here, I will use this to provide information about what is
% stored in each of the matrices in the xp.data cell array. (Alternatively,
% we could also make each of these matrices an xPlt object!)
meta = struct;
meta.datainfo(1:2) = nDDictAxis;
meta.datainfo(1).name = 'time(ms)';
meta.datainfo(1).values = time;
meta.datainfo(2).name = 'cells';
meta.datainfo(2).values = [];
xp.meta = meta;
clear meta



%% Validate & get some properties of the data
xp.checkDims;       % Makes sure all the dimensions match up (e.g. xp.axis 
                    % must have the same length as number of dimensions in
                    % xp.data, and the size of each dimension must match
                    % the number of labels
                    
xp.fixAxes;         % This attempts to automatically fix any dimension
                    % mismatches between xp.data and the axis labels. This
                    % particular command does nothing because there are no
                    % errors in xp

                    
% Make a "bad" xPlt class, containing errors.
disp(size(xp.data));                % The 4th dimension of xp.data is of size 8
disp(xp.axis(4).values);            % It's corresponding axis should have 8 labels
xp_bad = xp; 
xp_bad.axis(4).values={'test'};     % Reduce this to 1 (mismatch)

% Check errors in new class (this produces an error, so disabling it)
xp_bad.checkDims;

% Auto fix errors in labels
xp_fixed = xp_bad.fixAxes; 

% View new labels
xp_fixed.axis(4).values         % The original cell array had been replaced by 
                                % one of appropriate length
                                
% View summary of the classes. Tells the dimensionality of xp.data and also
% the number of labels in each axis.
xp.getaxisinfo
xp_fixed.getaxisinfo

%% Indexing data

% Indexing works just like with normal data. This creates a new xPlt object
% based on the original data, with the correct axis labels
xp4 = xp(:,:,1,8);                  % ## Update - This does the same thing as xp.subset([],[],1,8), which was the old way
                                    % of pulling data subsets. Note that [] and : are equivalent.
xp4.getaxisinfo

% Similarly, can index string axes using regular expressions
% Pull out sodium mechs only
xp5 = xp(:,:,1,'iNa*');
xp5.getaxisinfo

% Pull out synaptic state variables
xp5 = xp(:,:,1,'_s');
xp5.getaxisinfo

% Lastly, you can reference xp.data with the following shorthand
% (This is the same as xp.data(:,:,1,8). Regular expressions dont work in this mode)
mydata = xp{:,:,1,8};              
mydata2 = xp.data(:,:,1,8);
disp(isequal(mydata,mydata2));

clear mydata mydata2 xp4 xp5

%% Plot for 2D parameter sweep
% Tip: don't try to understand what recursivePlot is doing - instead, try
% putting break points in the various function handles to see how this
% command works.

xp4 = xp(:,:,'E','v');
recursivePlot(xp4,{@xp_subplot,@xp_matrix_basicplot},{[1,2]},{{0,0},{}});



%% Plot for 3D data

xp4 = xp(:,:,:,8);
xp4.getaxisinfo

recursivePlot(xp4,{@xp_subplot_grid3D,@xp_matrix_basicplot},{[3,1,2]},{{},{}});




%% Plot for 3D data (larger one)

xp4 = xp(:,:,1,:);
xp4 = xp4.squeeze;
xp4.getaxisinfo

% recursivePlot(xp4,{@xp_subplot,@xp_subplot,@xp_matrix_basicplot},{1:2,3},{{[],1},{1,1},{}});
% recursivePlot(xp4,{@xp_subplot_grid3D,@xp_subplot,@xp_matrix_basicplot},{1:2,3},{{},{0,1},{}});
% recursivePlot(xp4,{@xp_subplot_grid3D,@xp_matrix_basicplot},{[3,1,2]},{{},{}});
recursivePlot(xp4,{@xp_subplot_grid3D,@xp_matrix_basicplot},{[3,1,2]},{{},{}});



%% Use nested images to fit 3D data into a 2D subplot

xp4 = xp(:,:,:,8);
xp4.getaxisinfo

recursivePlot(xp4,{@xp_subplot_grid3D,@xp_subplot_grid3D,@xp_matrix_basicplot},{[1,2,4],3},{{},{0,1},{}});
% recursivePlot(xp4,{@xp_subplot_grid3D,@xp_subplot,@xp_matrix_basicplot},{[1,2,4],3},{{},{0,1},{}});

%% Combine two xPlt objects

xp3 = xp(2,:,'E','v');
xp3.getaxisinfo

xp4 = xp(:,3,'E','v');
xp4.getaxisinfo

xp5 = merge(xp3,xp4);

recursivePlot(xp5,{@xp_subplot,@xp_matrix_basicplot},{[1,2]},{{0,0},{}});


%% Test packDims
% Analogous to cell2mat.
clear xp2 xp3 xp4 xp5

% Start by taking a smaller subset of the original xp object.
% xp2 = xp.subset(2,2,[],[1,3,5:8]);      % Selection based on index locations
xp2 = xp(2,2,:,'(v|^i||ISYN$)');  % Same thing as above using regular expression. Selects everything except the _s terms. "^" - beginning with; "$" - ending with
xp2 = xp2.squeeze;
xp2.getaxisinfo;

% Note that xp2 is sparse (there are some empty cells)
disp(xp2.data);         % (E cells don't receive AMPA synapses, and I cells don't receive GABAA synapses)

% Now pack dimension two (columns) of xp2 into xp2.data.
src = 2;                    % Take 2nd dimension in xp2
dest = 3;                   % Pack into 3rd dimension in xp2.data matrix
xp3 = xp2.packDim(src,dest);


% Check dimensionality of xp3.data
disp(xp3.data)             % The dimension "variables", which was dimension 2
                           % in xp2, is now dimension 3 in xp3.data.
                           % Now xp3.data is time x cells x variables

% View axis of xp3
xp3.getaxisinfo;            % The dimension "variables" is now missing

% However, the information is stored in the nDDictAxis matrix_dim_3, a field of xp3.meta.
xp3.meta.matrix_dim_3.getaxisinfo

% And if dimension 3 of each cell in xp3.data is unpacked using unpackDim,
% xp3.meta.matrix_dim_3 will be used to provide axis info for the new
% xPlt object.
xp4 = xp3.unpackDim(dest, src);
xp4.getaxisinfo;

% Unless new axis info is provided, that is.
xp4 = xp3.unpackDim(dest, src, 'New_Axis_Names');
xp4.getaxisinfo;

%% This seems to be broken.
xp4 = xp3.unpackDim(dest, src, 'New_Axis_Names', {'One','Two','Three','Four','Five','Six'});
xp4.getaxisinfo;

%% So is this.
xp4 = xp3.unpackDim(dest, src, 'New_Axis_Names', []);
xp4.getaxisinfo;

%% And this.
xp4 = unpackDim(xp3, dest, src, 'New_Axis_Names', 1:6);
xp4.getaxisinfo;

%%
% Note some of this data is sparse!
temp1 = squeeze(xp3.data{1}(100,:,:));  % Pick out a random time point
temp2 = squeeze(xp3.data{2}(100,:,:));  % Pick out a random time point
figure; 
subplot(211); imagesc(temp1);

ylabel('Cells');
xlabel(xp2.axis(2).name); 
set(gca,'XTick',1:length(xp2.axis(2).values)); set(gca,'XTickLabels',strrep(xp2.axis(2).values,'_',' '));


subplot(212); imagesc(temp2);
ylabel('Cells');
xlabel(xp2.axis(2).name); 
set(gca,'XTick',1:length(xp2.axis(2).values)); set(gca,'XTickLabels',strrep(xp2.axis(2).values,'_',' '));


%% Use packDim to average across cells

xp2 = xp;
xp2 = xp(:,:,:,'v');  % Same thing as above using regular expression. Selects everything except the _s terms. "^" - beginning with; "$" - ending with
xp2 = xp2.squeeze;
%
% Average across all cells
xp2.data = cellfun(@(x) mean(x,2), xp2.data,'UniformOutput',0);

% % Convert xp2.data from a matrix into an xPlt object as well. This is
% % useful for keeping track of axis names. 
% mat_ax_names = {'Time','Cell Number'};
% mat_ax_values = {1:10001, []};
% 
% % xp2.data = Cell_2_nDDict(xp2.data,mat_ax_names,mat_ax_values);

% Pack E and I cells together
src=3;
dest=2;
xp3 = xp2.packDim(src,dest);

%%
% Plot 
recursivePlot(xp3,{@xp_subplot_grid3D,@xp_matrix_basicplot},{[1,2]},{{},{}});


%% Use packDim to average over synaptic currents
% Analogous to cell2mat
% See also plotting material by Hadley Wickham

% First, pull out synaptic current variables
xp2 = xp(:,:,:,'(ISYN$)');  % Same thing as above using regular expression. Selects everything except the _s terms. "^" - beginning with; "$" - ending with
xp2.getaxisinfo;

% Second, put this into matrix form, so we can average over them
xp3 = xp2.packDim(4,3);
disp(xp3.data)              % xp3.data is now 3D, with the 3rd dim denoting synaptic current
xp3 = xp3.squeeze;
xp3.getaxisinfo;

% Average across membrane currents
xp3.data = cellfun(@(x) nanmean(x,3), xp3.data,'UniformOutput',0);

% Plot 
recursivePlot(xp3,{@xp_subplot_grid3D,@xp_matrix_basicplot},{[3,1,2]},{{},{}});

%% Test mergeDims
% Analogous to Reshape.

% This command combines two (or more) dimensions into a single dimension.
xp2 = xp.mergeDims([3,4]);
xp2.getaxisinfo;

%% Use mergeDims to convert to Jason's DynaSim format
% Analogous to Reshape

% This combines the 2D "vary" sweep into a single dimension. It also
% combines all populations and variables into a single 1D list. Thus, Axis
% 1 is equivalent to Jason's structure array - data(1:9). Axis 4 is
% equivalent to the structure fields in Jason's DynaSim structure.
xp2 = xp.mergeDims([3,4]);
xp2 = xp2.mergeDims([1,2]);
xp2.getaxisinfo;

% Squeeze out the empty dimensions.
xp3 = squeeze(xp2);
xp3.getaxisinfo;



%% Load nDDict structure of images

% Import plot files
file = 'demo_sPING_3b';
data = ImportPlots(file);

% Load into DynaSim structure
[data_table,column_titles] = DataField2Table (data,'plot_files');

% Preview the contents of this table
previewTable(data_table,column_titles);

% The entries in the first column contain the paths to the figure files.
% There can be multiple figures associated with each simulation, which is
% why these are cell arrays of strings.
disp(data_table{1}{1})
disp(data_table{1}{2})

% Import the linear data into an xPlt object
xp = xPlt;
X = data_table{1}; axislabels = data_table(2:end);
xp = xp.importLinearData(X, axislabels{:});
xp = xp.importAxisNames(column_titles(2:end));

recursivePlot(xp,{@xp_subplot_grid3D,@xp_plotimage},{[1,2]},{{},{.5}});


%% To implement
% 
% 
% Implement the following:
% + Make DynaSimPlotExtract more general
% + Starting work on PlotData2 - any new requests?
