function zipfname = zipDemoData(study_dir)
%zipDemoData - Zips demo data (utility function)
%
% Purpose: Packages demo data into zip file and stores it in its
% appropriate location. This is useful if you want to create demo data to
% package with DynaSim.
%
% Usage:
%   zipfname = zipDemoData(study_dir)
%
% Inputs:
%   study_dir: input 
%
% Outputs:
%   zipfname: string with location of zipped file
%
% Examples:
%   zipfname = zipDemoData('/Users/davestanley/Documents/DynaSimData/demos/demo_sPING_100cells_3x3')
%
% See also: unZipDemoData, demos_generate_data.m (demo script)

    ds_path = getDsVar('ds_path');
    demos_path = getDsVar('demos_path');
    demo_zips_path = getDsVar('demos_zips_path');

    mkdirSilent(demo_zips_path);

    % Remove prefix
    [~,study_name] = fileparts(study_dir);


    zipfname = fullfile(demo_zips_path,[study_name '.zip']);
    zip(zipfname, fullfile(demos_path,study_name))

end
