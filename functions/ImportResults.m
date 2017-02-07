function results = ImportResults(studyinfo,func)
% results = ImportResults(studyinfo,func)
% inputs:
%   - studyinfo: DynaSim studyinfo structure or study_dir
%   - func: function handle of analysis function whose results to return
% outputs:
%   - results: structure of results [num_sims x num_calls]
% todo: This command breaks when "results" are figures e.g. outputs of
% PlotData (dave, Feb 2017). Does not know how to "load" an image, nor does
% it recognize the image extensions. I wrote "ImportPlots" as a way around
% this, but there might be better solutions for differentiating "plots" from
% other "results"

if ischar(studyinfo) && isdir(studyinfo) % study directory
  study_dir=studyinfo;
  clear studyinfo
  studyinfo.study_dir=study_dir;
end
if isstruct(studyinfo) && isfield(studyinfo,'study_dir')
  % retrieve most up-to-date studyinfo structure from studyinfo.mat file
  studyinfo=CheckStudyinfo(studyinfo.study_dir);
  % get list of data_files from studyinfo
  result_functions=studyinfo.simulations(1).result_functions;
  matches=cellfun(@(x)isequal(x,func),result_functions);
  result_files=cellfun(@(x)x(matches),{studyinfo.simulations.result_files},'uni',0);
  num_instances=cellfun(@length,result_files);
  num_sims=length(studyinfo.simulations);
  for s=1:num_sims
    for i=1:num_instances
      if exist(result_files{s}{i},'file')
        load(result_files{s}{i},'result');
        results(s,i)=result;
      end
    end
  end
end
