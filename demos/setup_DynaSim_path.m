
%% Attempts to set up DynaSim path from within demos folder
% Author David Stanley stanleyd@bu.edu

% If DynaSim isn't already in the system path, this code attempts to guess
% the location of the main DynaSim folder. If cannot find it, returns
% an error.
if ~exist('SimulateModel','file')                                                               % If DynaSim not in path...  
    [parentfolder_full,~] = fileparts(pwd); [~,parentfolder] = fileparts(parentfolder_full);    % Get name of parent folder (perhaps there is a better way to do this!)
    [currfolder_full] = pwd; [~,currfolder] = fileparts(currfolder_full);                       % Get name of current folder
    if ~isempty(strfind(parentfolder,'DynaSim'))                                                % Check if parent folder is DynaSim
        addpath(genpath(parentfolder_full));                                                    %     If is, add to path
        rmpath(genpath(fullfile(parentfolder_full,'.git')));                                    %        (Remove .git folder)
    elseif ~isempty(strfind(currfolder,'DynaSim'))                                              % Check if current folder is DynaSim
        addpath(genpath(currfolder_full));                                                      %     If is, add to path
        rmpath(genpath(fullfile(currfolder_full,'.git')));                                      %         (Remove .git folder)
    else
        % Otherwise, return error
        error('Add the DynaSim folder to the MATLAB path - e.g. run addpath(genpath(DynaSimPath))');
    end
end
