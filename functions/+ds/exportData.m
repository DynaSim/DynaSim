function exportData(data,varargin)
%EXPORTDATA - export DynaSim data structure in various formats.
%
% Usage:
%   ds.exportData(data,varargin)
%
% Inputs:
%   - data: DynaSim data structure (see ds.checkData)
%   - options:
%     'filename'    : name of output data file (default: 'data.mat')
%     'format'      : mat. todo: csv, HDF. (default: 'mat')
%     'verbose_flag': whether to print log info (default: 0)
% 
% See also: dsImport, ds.checkData, dsSimulate

options=ds.checkOptions(varargin,{...
  'filename','data.mat',[],... % name of output data file
  'format','mat',[],... % mat. todo: csv, HDF
  'verbose_flag',0,{0,1},... % whether to print log info
  },false);

switch lower(options.format)
  case 'mat'
    if numel(data)==1
      % split DynaSim data structure into separate variables saved to a
      % v7.3 mat-file for subsequent HDF-style loading with matfile()
      vars=fieldnames(data);
      for i=1:length(vars)
        eval(sprintf('%s=data.%s;',vars{i},vars{i}));
      end
      save(options.filename,vars{:},'-v7.3');
    else
      save(options.filename,'data','-v7.3');
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
