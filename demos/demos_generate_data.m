%% Generate data
% This is a script for generating small sets of dat for use in testing. It
% saves the data in .mat files, which seems to be a bit more efficient than
% saving as a DynaSim path


%% Set up paths 
% Get ready...

% Format
format compact
clear all
restoredefaultpath

% Check if in right folder
[parentfolder,currfolder] = fileparts(pwd);
if ~strcmp(currfolder,'demos'); error('Should be in demos folder to run this code.'); end

% add DynaSim toolbox to Matlab path
addpath(genpath(parentfolder)); % comment this out if already in path

% Study directory
output_directory = fullfile(parentfolder,'outputs');
study_dir = fullfile(output_directory,'demo_sPING_3b_2plots');


%% Generate a huge dataset
% Run simulation - Sparse Pyramidal-Interneuron-Network-Gamma (sPING)

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
  'E'   ,'Iapp',[0:10:10];      % amplitude of tonic input to E-cells
  %'I'   ,'Iapp',[0 5 10];      % amplitude of tonic input to E-cells
  'I->E','tauD',[5:5:5]       % inhibition decay time constant from I to E
  };
SimulateModel(s,'save_data_flag',1,'study_dir',study_dir,...
                'vary',vary,'verbose_flag',1, 'downsample_factor', 10, ...
                'save_results_flag',1,'plot_functions',{@PlotData,@PlotData},'plot_options',{{'format','png','visible','off','figwidth',0.5,'figheight',0.5}, ...
                {'format','png','visible','off','plot_type','rastergram','figwidth',0.5,'figheight',0.5}} );

%% Load the data

% Load data in traditional DynaSim format
data=ImportData(study_dir);

% Import the data images

% Import plot files
data_img = ImportPlots(study_dir);


%% Resize data as needed

Num_cells_to_keep = 20;
downsample_factor = 2;
mydata = xp.data;

for i = 1:numel(data)
    if ~isempty(mydata{i})
        %mydata{i} = mydata{i}(:,1:Num_cells_to_keep);
        mydata{i} = mydata{i}(:,1:round(size(mydata{i},2)/downsample_factor));
    end
end

xp.data = mydata;

data = xPlt2DynaSim(xp);


save('sample_data_dynasim_2plots.mat','data','data_img');

%% 


% Load into DynaSim structure
[data_table,column_titles] = DataField2Table (data_img,'plot_files');

% Preview the contents of this table
previewTable(data_table,column_titles);

% The entries in the first column contain the paths to the figure files.
% There can be multiple figures associated with each simulation, which is
% why these are cell arrays of strings.
disp(data_table{1}{1})
disp(data_table{1}{2})

% Import the linear data into an xPlt object
xp_img = xPlt;
X = data_table{1}; axislabels = data_table(2:end);
xp_img = xp_img.importLinearData(X, axislabels{:});
xp_img = xp_img.importAxisNames(column_titles(2:end));



