%% % % % % % % % % % % % % % % SETTING UP xPlt OBJECT % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%% Set up paths 
% Get ready...

% Format
format compact

% Check if in right folder
[parentfolder,currfolder] = fileparts(pwd);
if ~strcmp(currfolder,'demos'); error('Should be in demos folder to run this code.'); end

% Add DynaSim to path if it's not already there
if exist('setup_DynaSim_path','file')
    setup_DynaSim_path;
else
    error('Add the DynaSim folder to the MATLAB path - e.g. run addpath(genpath(DynaSimPath))');
end

fprintf('Note1 - I am planning to rename nDDict to MDD (MultiDimensional Dictionary).\n');
fprintf('Note2 - I have moved the rest of this demos script to the ../MDD now in\n ');
fprintf('order for it to updated in sync with the MDD repo.\n');

% Open demos_xPlt.m script
edit ../dependencies/MDD/demos_xPlt.m
cd  ../dependencies/MDD/

