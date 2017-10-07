function results = dsImportResults(studyinfo,func, varargin)
%IMPORTRESULTS - Import analysis result of a simulation
%
% Usage:
%   results = dsImportResults(studyinfo,func)
%
% Inputs:
%   - studyinfo: DynaSim studyinfo structure or study_dir
%   - func: function handle of analysis function whose results to return
%
% Outputs:
%   - results: structure of results [num_sims x num_calls]
%
% TODO:
%   - This command breaks when "results" are figures e.g. outputs of dsPlot
%   (dave, Feb 2017). Does not know how to "load" an image, nor does it
%   recognize the image extensions. I wrote "dsImportPlots" as a way around this,
%   but there might be better solutions for differentiating "plots" from other
%   "results"
% 
% Author: Jason Sherfey, PhD <jssherfey@gmail.com>
% Copyright (C) 2016 Jason Sherfey, Boston University, USA


if ischar(studyinfo) && isdir(studyinfo) % study directory
  study_dir=studyinfo;
  clear studyinfo
  studyinfo.study_dir=study_dir;
end

if isstruct(studyinfo) && isfield(studyinfo,'study_dir')
  % retrieve most up-to-date studyinfo structure from studyinfo.mat file
  studyinfo=dsCheckStudyinfo(studyinfo.study_dir, varargin{:});
  if exist('study_dir','var')
    studyinfo.study_dir=study_dir;
  end
  
  % get list of data_files from studyinfo
  result_functions=studyinfo.simulations(1).result_functions;
  matches=cellfun(@(x) strcmp(func2str(x), func2str(func)),result_functions);
  if ~any(matches)
    wprintf('Didnt find match for result function parameter')
    return
  end
  result_files=cellfun(@(x)x(matches),{studyinfo.simulations.result_files},'uni',0);
  num_instances=cellfun(@length,result_files);
  num_sims=length(studyinfo.simulations);
  
%   result_files_exist=cellfun(@(x) cellfun(@exist, x),result_files)==2;
%   if ~any(result_files_exist)
%     % convert original absolute paths to paths relative to study_dir
%     for s=1:num_sims
%       for i=1:num_instances(s)
%         [~,fname,fext]=fileparts2(result_files{s}{i});
%         result_files{s}{i}=fullfile(studyinfo.study_dir,'data',[fname fext]);
%       end
%     end
%   end
  
  for s=1:num_sims
    for i=1:num_instances(s)
      %check absolute path
      if exist(result_files{s}{i},'file')
        load(result_files{s}{i},'result');
        results(s,i)=result;
      else
        %check relative path
        [~,fname,fext]=fileparts2(result_files{s}{i});
        result_files{s}{i}=fullfile(studyinfo.study_dir,'data',[fname fext]);
        if exist(result_files{s}{i},'file')
          load(result_files{s}{i},'result');
          results(s,i)=result;
        end
      end
    end
  end
end
