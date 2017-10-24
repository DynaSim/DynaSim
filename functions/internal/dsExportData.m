function dsExportData(data,varargin)
%EXPORTDATA - export DynaSim data structure in various formats.
%
% Usage:
%   dsExportData(data,varargin)
%
% Inputs:
%   - data: DynaSim data structure (see dsCheckData)
%   - options:
%     'filename'    : name of output data file (default: 'data.mat')
%     'format'      : mat. todo: csv, HDF. (default: 'mat')
%     'verbose_flag': whether to print log info (default: 0)
%
% See also: dsImport, dsCheckData, dsSimulate
%
% Author: Jason Sherfey, PhD <jssherfey@gmail.com>
% Copyright (C) 2016 Jason Sherfey, Boston University, USA

options=dsCheckOptions(varargin,{...
  'filename','data.mat',[],... % name of output data file
  'format','mat',[],... % mat. todo: csv, HDF
  'matCompatibility_flag',1,{0,1},... % whether to save mat files in compatible mode, or to prioritize > 2GB VARs
  'verbose_flag',0,{0,1},... % whether to print log info
  },false);

switch lower(options.format)
  case 'mat'
    if numel(data)==1
      % split DynaSim data structure into separate variables saved to a
      % mat-file for subsequent loading with matfile()
      vars=fieldnames(data);
      for i=1:length(vars)
        eval(sprintf('%s=data.%s;',vars{i},vars{i}));
      end
      if options.matCompatibility_flag
        try
          save(options.filename,vars{:},'-v7');
        catch
          fprintf('Data is not ''-v7'' compatible. Setting ''matCompatibility_flag'' to 0.\n')
          options.matCompatibility_flag = 0;
          if strcmp(reportUI,'matlab')
            save(options.filename,vars{:},'-v7.3');
          else
            save(options.filename,vars{:},'-hdf5'); % hdf5 format in Octave
          end
        end
      else
        if strcmp(reportUI,'matlab')
          save(options.filename,vars{:},'-v7.3');
        else
          save(options.filename,vars{:},'-hdf5'); % hdf5 format in Octave
        end
      end
    else
      if options.matCompatibility_flag
        try
          save(options.filename,'data','-v7');
        catch
          fprintf('Data is not ''-v7'' compatible. Setting ''matCompatibility_flag'' to 0.\n')
          options.matCompatibility_flag = 0;
          if strcmp(reportUI,'matlab')
            save(options.filename,'data','-v7.3');
          else
            save(options.filename,'data','-hdf5'); % hdf5 format in Octave
          end
        end
      else
        if strcmp(reportUI,'matlab')
          save(options.filename,'data','-v7.3');
        else
          save(options.filename,'data','-hdf5'); % hdf5 format in Octave
        end
      end
    end
  case 'csv'
    % write csv file
    % ExportCSV(data,options.filename);
  case 'hdf'
    % ExportHDF(data,options.filename);
  otherwise
    % not recognized
end

if options.verbose_flag
  fprintf('\tData saved to %s\n',options.filename);
end
