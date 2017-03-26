machine='dave';

% % Get home folder
if ispc;
    userdir= getenv('USERPROFILE'); 
    docsfolder = 'Documents';           % Is this correct on PC?
    warning('Double check docs folder on PC is correct');
else
    userdir= getenv('HOME');
    docsfolder = 'Documents';
end

% % System-specific customization
demos_output = fullfile(userdir,docsfolder,'DynaSimData','demos');

