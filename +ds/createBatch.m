function studyinfo = createBatch(base_model,modifications_set,varargin)
%CREATEBATCH - create and submit jobs to run sets of simulations or analyses.
%
% Usage:
%   studyinfo=ds.createBatch(model,varargin)
%
% Inputs:
%   - DynaSim model (see ds.generateModel)
%   - modifications_set (as returned by ds.vary2Modifications())
%     - options:
%       'compile_flag'  : whether to compile simulation using coder instead of
%                         interpreting Matlab {0 or 1} (default: 0)
%       'verbose_flag'  : whether to display informative messages/logs (default: 0)
%       'overwrite_flag': whether to overwrite existing data files {0 or 1} (default: 0)
%     - options for cluster computing:
%       'sims_per_job'  : number of simulations to run per cluster job (default: 1)
%       'memory_limit'  : memory to allocate per cluster job (default: '8G')
%       'batch_dir'     : where to save job scripts
%       'qsub_mode'     : whether to use SGE -t array for 1 qsub, mode: 'array'; or
%                         qsub in csh for loop, mode: 'loop'. (default: 'loop').
%     - options for parallel computing: (requires Parallel Computing Toolbox)
%       - Note: parallel computing has been DISABLED for debugging...
%       'parallel_flag' : whether to use parfor to run simulations {0 or 1} (default: 0)
%       'num_cores'     : number of cores to specify in the parallel pool
%
% Outputs:
%   - DynaSim studyinfo (see ds.checkStudyinfo for schema info)
%
% Dependencies: ds.setupStudy, ds.updateStudy
%
% See also: ds.generateModel, ds.simulateModel, ds.checkStudyinfo, ds.vary2Modifications

% check inputs
options=ds.checkOptions(varargin,{...
  'compile_flag',0,{0,1},...
  'parallel_flag',0,{0,1},...
  'num_cores',4,[],... % # cores for parallel processing (SCC supports 1-12)
  'sims_per_job',1,[],... % how many sims to run per cluster job
  'memory_limit','8G',[],... % how much memory to allocate per cluster job
  'batch_dir',[],[],...
  'simulator_options',[],[],...
  'verbose_flag',0,{0,1},...
  'overwrite_flag',0,{0,1},...
  'process_id',[],[],... % process identifier for loading studyinfo if necessary
  'qsub_mode','loop',{'loop','array'},... % whether to submit jobs as an array using qsub -t or in a for loop
  'one_solve_file_flag',0,{0,1},... % use only 1 solve file of each type, but can't vary mechs yet
  'solver','rk4',{'euler','rk1','rk2','rk4','modified_euler','rungekutta','rk','ode23','ode45',...
    'ode1','ode2','ode3','ode4','ode5','ode8','ode113','ode15s','ode23s','ode23t','ode23tb'},... % DynaSim and built-in Matlab solvers
  },false);

% Set up studyinfo structure, study directory and output file names
%[studyinfo,options.simulator_options]=ds.setupStudy(base_model,modifications_set,options.simulator_options);
[studyinfo,options.simulator_options]=ds.setupStudy(base_model,'modifications_set',modifications_set,'simulator_options',options.simulator_options,'process_id',options.process_id);
study_file=fullfile(studyinfo.study_dir,'studyinfo.mat');
num_simulations=length(modifications_set);

% check whether study has already completed
% if options.overwrite_flag==0
  [~,s]=ds.monitorStudy(studyinfo.study_dir,'verbose_flag',0,'process_id',options.process_id);
  if s==1 % study finished
    if options.verbose_flag
      fprintf('Study already finished. Not creating new batch.\n');
    end
    
    studyinfo=ds.checkStudyinfo(studyinfo.study_dir,'process_id',options.process_id);
    
    return;
  end
% end

% Check if attempting to compile Experiment
if isa(options.simulator_options.experiment,'function_handle') && options.compile_flag
  options.compile_flag=0;
  options.simulator_options.compile_flag=0;
  wprintf('solver compilation is not available for experiments run on the cluster.');
end

% create base solve_file (m- or MEX-file)
solve_file = ds.getSolveFile(base_model,studyinfo,options.simulator_options);

% add appropriate MEX extension if compiled
%   NOTE: the extension depends on whether a 32-bit or 64-bit machine is used
if options.compile_flag
  % get directory where solve_file is located and the file name
  [fpath,fname]=fileparts(solve_file);
  
  % get list of files in solve directory
  D=dir(fpath);
  
  % get extension of files with matching names
  tmp=regexp({D.name},[fname '\.(\w+)$'],'tokens','once');
  tmp=unique([tmp{:}]);
  
  % append extension to solve_file if regular mex (not non supported matlab solver mex)
  if ~any(strcmp(options.solver, {'ode113','ode15s','ode23s','ode23t','ode23tb'})) % not mex supported)
    full_solve_file=[solve_file '.' tmp{1}];
  
    % remove '_mex' suffix from solve_file for compatibility with ds.getSolveFile()
    solve_file=regexp(solve_file,'(.+)_mex$','tokens','once');
    solve_file=solve_file{1};
  else
    full_solve_file=solve_file;
  end
else
  full_solve_file=solve_file;
end
studyinfo.base_solve_file=solve_file;

% set name of batch_dir for this study
% get study-specific timestamp
timestamp = datestr(studyinfo.study_id,'yyyymmddHHMMSS');

% get home directory of the current user
[o,home]=system('echo $HOME');

% create batch directory
[study_dir_path,study_dir_name,study_dir_suffix]=fileparts(studyinfo.study_dir);
study_dir_name=[study_dir_name study_dir_suffix];
batch_dir = fullfile(strtrim(home),'batchdirs',study_dir_name);
%batch_dir = fullfile(strtrim(home),'batchdirs',['Batch' timestamp]);

if options.verbose_flag
  fprintf('\nPREPARING BATCH:\n');
end

% create batch_dir to hold m-file jobs and more for this study
if ~exist(batch_dir,'dir')
  if options.verbose_flag
    fprintf('Creating batch jobs directory: %s\n',batch_dir);
  end
  mkdir(batch_dir);
end

% copy study_file to batchdir
[success,msg]=copyfile(study_file,batch_dir);
if ~success
  error(msg)
end

% locate DynaSim toolbox
dynasim_path=fileparts(fileparts(which(mfilename))); % root is one level up from directory containing this function
dynasim_functions=fullfile(dynasim_path,'functions');

% locate mechanism files
[mech_paths,mech_files]=ds.locateModelFiles(base_model);

% add paths to studyinfo structure
studyinfo.paths.dynasim_functions = dynasim_path;
studyinfo.paths.mechanisms = mech_paths;
studyinfo.paths.batch_dir = batch_dir;

% edit simulator_options so that subsequent single simulations do not create further batches
options.simulator_options.parallel_flag = 0;
options.simulator_options.cluster_flag = 0;
options.simulator_options.verbose_flag = 1; % turn on verbose for cluster logs (stdout)

% create jobs (k)
jobs={};
skip_sims=[];
sim_cnt=0;
num_jobs=ceil(num_simulations/options.sims_per_job);

if options.verbose_flag && ~options.one_solve_file_flag
  fprintf('Creating %g cluster jobs in batch directory for %g simulations...\n',num_jobs,num_simulations);
end

if ~options.one_solve_file_flag
  for k = 1:num_jobs
    job_file=fullfile(batch_dir,sprintf('sim_job%g.m',k));
    sim_ids=sim_cnt+(1:options.sims_per_job);
    sim_ids=sim_ids(sim_ids<=num_simulations);
    sim_cnt=sim_cnt+options.sims_per_job;
    
    % only run simulations if data_file does not exist (unless overwrite_flag=1)
    proc_sims=[]; % simulations to run in this job
    for j=1:length(sim_ids)
      %if ~exist(studyinfo.simulations(sim_ids(j)).data_file,'file') || options.overwrite_flag
      if (options.simulator_options.save_data_flag && ~exist(studyinfo.simulations(sim_ids(j)).data_file,'file')) || ...
          (options.simulator_options.save_results_flag && ~exist(studyinfo.simulations(sim_ids(j)).result_files{1},'file')) || ...
          options.overwrite_flag
        proc_sims=[proc_sims sim_ids(j)];
      end
    end
    
    if isempty(proc_sims)
      skip_sims=[skip_sims sim_ids];
      continue; % skip this job
    else
      skip_sims=[skip_sims setdiff(sim_ids,proc_sims)];
    end
    
    WriteSimJob(proc_sims,job_file); % create job script
    
    jobs{end+1}=job_file;
    
    % add job_file to studyinfo
    for j=1:length(proc_sims)
      studyinfo.simulations(proc_sims(j)).job_file=job_file;
    end
      
  end %for
else %one_solve_file_flag
  job_file=fullfile(batch_dir, 'sim_job.m'); %always use this file
  
  WriteSimJob([], job_file); % create job script
  
  % TODO: compile job_file
%   if options.compile_flag
%     job_file = ds.prepareMEX(job_file, options); %compile the job file
%   end
  
  % add job_file to studyinfo
  [studyinfo.simulations.job_file] = deal(job_file); %same job file for all sims
end %if

% TODO: have job scripts create lock for each sim; then, check for lock
%   file in this function and skip simulations if they're already running
%   (similar to skipping if data_file already exists).

% exit if there are no simulations to run
if numel(skip_sims)==num_simulations
  if options.verbose_flag
    fprintf('Data exists for all simulations in study. nothing to do.\n');
  end
  
  return;
end

% create script_filename (list of vars or jobs)
if strcmp(options.qsub_mode, 'loop')
  script_filename = fullfile(batch_dir,'scriptlist.txt');
  if options.verbose_flag
    fprintf('Creating file listing jobs in batch directory: %s\n',script_filename);
  end
end

% write to file
if strcmp(options.qsub_mode, 'loop')
  fScript=fopen(script_filename,'wt');
  
  for j=1:length(jobs)
    [~,thisFilename]=fileparts(jobs{j});
    fprintf(fScript,'%s\n',thisFilename);
  end
  
  fclose(fScript);
end

if ~options.one_solve_file_flag
  % copy solve_file to each sim-specific solve sub-directory: <study_dir>/solve/sim<k>/<solve_file>
  %   and set sim-specific studyinfo
  if options.verbose_flag
    fprintf('Creating distinct solver sub-directories for %g simulations...\n',num_simulations);
  end
  %[solve_path,solve_name,solve_ext]=fileparts(solve_file);
  [solve_path,solve_name,solve_ext]=fileparts(full_solve_file);

  % prepare studyinfo with simulation-specific metadata
  for sim=1:num_simulations
    if ismember(sim,skip_sims)
      continue; % do nothing with this simulation
    end

    this_solve_path=fullfile(solve_path,['sim' num2str(sim)]);

    if ~exist(this_solve_path,'dir')
      %if options.verbose_flag
      %  fprintf('creating solver sub-directory for simulation #%g: %s\n',sim,this_solve_path);
      %end
      mkdir(this_solve_path);
    end

    [success,msg]=copyfile(full_solve_file,this_solve_path); % creates solve sub-directory and copies the base solve file to it

    if ~success, error(msg); end

    % set studyinfo solve_file to use for this simulation
    this_solve_file=fullfile(this_solve_path,[solve_name solve_ext]);
    studyinfo.simulations(sim).solve_file=this_solve_file;
    studyinfo.simulations(sim).batch_dir=batch_dir;
    studyinfo.simulations(sim).simulator_options=options.simulator_options;

    % set solve_file for this simulation (will be passed to ds.simulateModel as an option)
    if isa(options.simulator_options.experiment,'function_handle')
      studyinfo.simulations(sim).simulator_options.solve_file=[];
    else
      studyinfo.simulations(sim).simulator_options.solve_file=this_solve_file;
    end

    % set vary to [] to avoid each sim expanding to a new set
    studyinfo.simulations(sim).simulator_options.vary=[];
    studyinfo.simulations(sim).simulator_options.sim_id=studyinfo.simulations(sim).sim_id;
    studyinfo.simulations(sim).error_log='';
    
    % copy studyinfo to batch_dir for each simulation
    %this_study_file=fullfile(batch_dir,sprintf('studyinfo_%g.mat',sim));
    %save(this_study_file,'studyinfo');
    
    % copy studyinfo file to batch_dir for each simulation
    for sim=1:num_simulations
      this_study_file=fullfile(batch_dir,sprintf('studyinfo_%g.mat',sim));
      
      if sim==1
        save(this_study_file,'studyinfo');
        first_study_file=this_study_file;
      else
        % use copyfile() after saving first b/c >10x faster than save()
        [success,msg]=copyfile(first_study_file,this_study_file);
        
        if ~success, error(msg); end
      end
    end

  end %sim
else %one_solve_file_flag
  % set studyinfo solve_file to use for this simulation
  [studyinfo.simulations.solve_file] = deal(full_solve_file);
  [studyinfo.simulations.batch_dir] = deal(batch_dir);
  [studyinfo.simulations.simulator_options] = deal(options.simulator_options);

  % set vary to [] to avoid each sim expanding to a new set
  for iSim = 1:num_simulations
    studyinfo.simulations(iSim).simulator_options.vary = [];
    studyinfo.simulations(iSim).simulator_options.sim_id = studyinfo.simulations(iSim).sim_id;
  end
  [studyinfo.simulations.error_log] = deal('');
 
  % copy studyinfo file to batch_dir since more information now
  batch_study_file = fullfile(batch_dir,'studyinfo.mat');
  save(batch_study_file,'studyinfo');
end

% update studyinfo on disk
ds.studyinfoIO(studyinfo,study_file,options.simulator_options.sim_id,options.verbose_flag);

% TODO: remove old lock files ..

% check for qsub on system
[status,result]=system('which qsub');
if isempty(result)
  [~,host] = system('hostname');
  fprintf('qsub not found on host (%s).\n',strtrim(host));
  fprintf('Jobs NOT submitted to cluster queue.\n');
else % on cluster with qsub
  % submit jobs (e.g., fScripjobs_memlimit batch_dir 32G)
  % check status of study
  [~,s]=ds.monitorStudy(studyinfo.study_dir,'verbose_flag',0);
  
  if s~=1 % study not finished
    % submit jobs using shell script
    [batch_dir_path,batch_dir_name,batch_suffix]=fileparts(batch_dir);
    batch_dir_name=[batch_dir_name batch_suffix];
    dsFnPath = fileparts(mfilename('fullpath'));
    
    if options.parallel_flag
      cmd=sprintf('qmatjobs_pct %s %s %g',batch_dir_name,options.memory_limit,options.num_cores);
      %status=system('which qmatjobs_pct');
    else
      if strcmp(options.qsub_mode, 'array') && ~options.one_solve_file_flag
        % TODO: remove old error and output files; put e and o in their own dirs
        cmd = sprintf('echo ''%s/qmatjob_array %s sim_job'' | qsub -V -hard -l ''h_vmem=%s'' -wd %s -N %s_sim_job -t 1-%i',...
          dsFnPath, batch_dir, options.memory_limit, batch_dir, batch_dir_name, num_jobs);
      elseif strcmp(options.qsub_mode, 'array') && options.one_solve_file_flag
        [~, job_filename] = fileparts(job_file); %remove path and extension
        cmd = sprintf('echo ''%s/qmatjob_array_one_file %s %s'' | qsub -V -hard -l ''h_vmem=%s'' -wd %s -N %s_sim_job -t 1-%i:%i',...
          dsFnPath, batch_dir, job_filename, options.memory_limit, batch_dir, batch_dir_name, num_simulations, options.sims_per_job);
        % NOTE: using num_simulations, not num_jobs, since the job_file will
        %   determine it's own sims to run
      elseif strcmp(options.qsub_mode, 'loop')
        cmd = sprintf('qmatjobs_memlimit_loop %s %s',batch_dir_name,options.memory_limit);
      end
      %status=system('which qmatjobs_memlimit');
    end
    
    % add shell script to linux path if not already there
    %if status~=0
      setenv('PATH', [getenv('PATH') ':' dynasim_functions]);
    %end
    
    if options.verbose_flag
      fprintf('Submitting cluster jobs with shell command: "%s"\n',cmd);
    end
    
    [status,result]=system(cmd);
    if status > 0
      fprintf('Submit command failed: "%s"\n',cmd);
      disp(result);
      return;
    else
      fprintf('Submit command status: \n');
      disp(result);
    end
    
    %if options.verbose_flag
      %fprintf('%g jobs successfully submitted!\ntip: use ds.monitorStudy(''%s'') or ds.monitorStudy(studyinfo) to track status of cluster jobs.\n',num_jobs,studyinfo.study_dir);
      fprintf('%g jobs successfully submitted!\n',num_jobs);
    %end
  elseif s==1 % study is finished
    if options.verbose_flag
      fprintf('Study already finished. not submitting jobs.\n');
    end
  end
end

%% NESTED FUNCTIONS
  function WriteSimJob(sim_ids,job_file)
    % purpose: write m-file to job_file to run simulations of sim_ids
    
    % create job file
    fjob=fopen(job_file,'wt');
    
    % load studyinfo using helper function to avoid busy file errors
    %fprintf(fjob,'studyinfo=ds.checkStudyinfo(''%s'',''process_id'',%g);\n',study_file,sim_ids(1));
    %fprintf(fjob,'load(''%s'',''studyinfo'');\n',study_file);
    
    if ~options.one_solve_file_flag
      % function declaration
      fprintf(fjob, 'function %s\n\n', job_file);
      
      % set IDs of simulations to run
      fprintf(fjob,'SimIDs=[%s]\n',num2str(sim_ids));
    else %only 1 file
      % function declaration
      [~, job_filename] = fileparts(job_file); %remove path and extension
      fprintf(fjob, 'function %s(simIDstart, simIDstep, simIDlast)\n\n', job_filename);
      
      if options.compile_flag
        fprintf(fjob, 'assert(isa(simIDstart, ''double''));\n');
        fprintf(fjob, 'assert(isa(simIDstep, ''double''));\n');
        fprintf(fjob, 'assert(isa(simIDlast, ''double''));\n');
      end
      
      % set IDs of simulations to run
      fprintf(fjob,'SimIDs = simIDstart:min(simIDlast,simIDstart+simIDstep-1);\n');
    end
    
    % loop over and run simulations in this job
    if options.parallel_flag
      % set parallel computing options
      %fprintf(fjob,'started=0;\npool=gcp(''nocreate'');\n');
      %fprintf(fjob,'if isempty(pool), parpool(%g); started=1; end\n',options.num_cores);
      fprintf(fjob,'c=parcluster;\nsaveAsProfile(c,''local_%g'');\nparpool(c,%g);\n',k,options.num_cores);
      
      % use parfor loop
      fprintf(fjob,'parfor s=1:length(SimIDs)\n');
    else
      % use for loop
      fprintf(fjob,'for s=1:length(SimIDs)\n');
    end
    
    fprintf(fjob,'\tSimID=SimIDs(s);\n');
    
    % each job should have try-catch to capture studyinfo.simulations(k).error_log
    fprintf(fjob,'\ttry\n');
    
    % load studyinfo for this simulation
    if ~options.one_solve_file_flag
      fprintf(fjob,'\t\tload(fullfile(''%s'',sprintf(''studyinfo_%%g.mat'',SimID)),''studyinfo'');\n',batch_dir);
    else
      fprintf(fjob,'\t\tload(fullfile(''%s'',''studyinfo.mat''),''studyinfo'');\n',batch_dir);
    end
    % compare paths between compute machine and studyinfo startup
    fprintf(fjob,'\t\t[valid,message]=ds.checkHostPaths(studyinfo);\n');
    
    if ~options.one_solve_file_flag
      fprintf(fjob,'\t\tif ~valid\n\t\t  lasterr(message);\n\t\t  for s=1:length(SimIDs), ds.updateStudy(studyinfo.study_dir,''sim_id'',SimIDs(s),''status'',''failed''); end\n\t\t  continue;\n\t\tend\n');
    end
    
    % simulate model with proper modifications and options
    fprintf(fjob,'\t\tsiminfo=studyinfo.simulations(SimID);\n');
    fprintf(fjob,'\t\toptions=rmfield(siminfo.simulator_options,{''modifications'',''studyinfo'',''analysis_functions'',''plot_functions'',''sim_id''});\n');
    fprintf(fjob,'\t\tkeyvals=ds.options2Keyval(options);\n');
    fprintf(fjob,'\t\tfprintf(''-----------------------------------------------------\\n'');\n');
    fprintf(fjob,'\t\tfprintf(''Processing simulation %%g (%%g of %%g in this job)...\\n'',SimID,s,length(SimIDs));\n');
    fprintf(fjob,'\t\tfprintf(''-----------------------------------------------------\\n'');\n');
    fprintf(fjob,'\t\tdata=ds.simulateModel(studyinfo.base_model,''modifications'',siminfo.modifications,''studyinfo'',studyinfo,''sim_id'',SimID,keyvals{:});\n');
    fprintf(fjob,'\t\tfor i=1:length(siminfo.result_functions)\n');
    fprintf(fjob,'\t\t\tds.analyzeData(data,siminfo.result_functions{i},''result_file'',siminfo.result_files{i},''save_data_flag'',1,siminfo.result_options{i}{:});\n');
    fprintf(fjob,'\t\tend\n');
    
    % add error handling
    fprintf(fjob,'\tcatch err\n');
    
    %fprintf(fjob,'\t\ttry delete(fullfile(studyinfo.study_dir,[''.locked_'' num2str(SimID)])); end\n');
    fprintf(fjob,'\t\tdisplayError(err);\n');
    fprintf(fjob,'\tend\n');
    
    % end loop over simulations in this job
    fprintf(fjob,'end\n');
    if options.parallel_flag
      %fprintf(fjob,'if started, delete(gcp); end\n');
      fprintf(fjob,'delete(gcp)\n');
    end
    
    % exit script
    [status,result]=system('which qsub');
    if ~isempty(result) % on cluster
      fprintf(fjob,'exit\n');
    end
    
    % close job file
    fclose(fjob);
  end % WriteSimJob

end

% -------------------
% in ds.createBatch():
% -------------------
% create solve_file
% for each sim k
%   copy to this_solve_file = /solve/sim<k>/solve_file
%   set studyinfo.simulations(k).solve_file=this_solve_file
% ...
% in sim job.m:
% ds.simulateModel(model{k},'solve_file',solve_file{k},...)
% where solve_file{k}=studyinfo.simulations(k).solve_file
