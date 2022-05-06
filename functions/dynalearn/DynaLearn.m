classdef DynaLearn < matlab.mixin.SetGet

    %
    % DynaModel; variational version
    % An integrated class for interactive and trainable DynaSim network models
    % Current learning method: Reinforcement learning - Delta rule
    %
    
    properties
        
        dlModel = []; % Dynasim's model struct, containing model's populations (layers) and connections, mechanisms ... .
        dlStudyDir = 'models/'; % Output data of last simulation.
        dlErrorsLog = []; % Log of errors since first trial.
        dlConnections = []; % List of connection names.
    
        dlTrialNumber = 0; % Last trial; how many times this model have been trianed with stimuli.
        dlLastOutputs = []; % Last output (eg. spike vector) of this model.
        dlLastError = 0; % Last error of this model for target, equivalent to the last value of errors_log property
        dlOutputLog = []; % Last expected output that this model should've generated.
        
        dsData = []; % Last simulation outputs
        dlOutputs = []; % Mex outputs
        dlVariables = []; % Mex variable labels
        dlMexFuncName = []; % Name of Mex function (e.g **********_mex.mex64
        
        dlPath = []; % Path which contains params.mat, mexfuncs, solver ...
        dlPathToFile = 'DynaSim/models/dlBaseModel';
        dlBaseVoltage = -77.4;
        dldT = .01; % Time step in ODEs (dt)
        
        dlDownSampleFactor = 10; 
        dlOptimalError = 1e9;
        dlUpdateError = 0;
        dlLastLambda = 1e-3;
        
        dlDeltaRatio = 1;
        dlLastDelta = 1;
        
    end
    
    methods

        function obj = DynaLearn(varargin) % Constructors, will be expanded
            
            disp("Creating Dyna model object ... ");
            set(obj, 'dlPathToFile', 'DynaSim/models/dlBaseModel');
            
            if nargin == 0
                
                obj = obj.dlLoad(obj.dlPathToFile);
                
            elseif nargin == 1

                model_ = varargin{1};
                set(obj, 'dlModel', model_);
                obj.dlInit(obj.dlStudyDir);
                set(obj, 'dlConnections', obj.dlGetConnectionsList());
                    
            elseif nargin == 2
                    
                model_ = varargin{1};  
                data_ = varargin{2};
%                 mexn_ = varargin{3};
                
                set(obj, 'dlModel', model_);
                set(obj, 'dlStudyDir', data_);
%                 set(obj, 'mex_func_name', mexn_);
                
                obj.dlInit(obj.dlStudyDir);
                set(obj, 'dlConnections', obj.dlGetConnectionsList());
                    
            else
                disp('Invalid use of DynaNet; pass a DynaSim struct and then address of parameters dataset file.');
            end
          
            [out, vars] = dsGetOutputList(obj.dlModel);
            set(obj, 'dlOutputs', out);
            set(obj, 'dlVariables', vars);
            obj.dlMexBridgeInit();
            
            disp("DynaLearn model created.");
            
        end

        function set.dlModel(obj, val) % Getter/setters
            
             if ~isstruct(val) 
                error('Model must be a struct');
             end
             obj.dlModel = val;
             
        end
        
        function set.dlStudyDir(obj, val)
            
             if ~strcmpi(class(val), 'string') && ~ischar(val)
                error('Study directory must be an address, a string');
             end
             obj.dlStudyDir = val;
             
        end
        
        function set.dlErrorsLog(obj, val)
            
             if ~strcmpi(class(val), 'double') 
                error('Errors log must be a double array');
             end
             obj.dlErrorsLog = val;
             
        end
        
        function set.dlConnections(obj, val)
            
             if ~iscell(val) 
                error('Connections must be a cell');
             end
             obj.dlConnections = val;
             
        end
        
        function set.dlTrialNumber(obj, val)
             
             obj.dlTrialNumber = floor(double(val));
             
        end
        
        function set.dlLastOutputs(obj, val)
             
             obj.dlLastOutputs = val;
             
        end
        
        function set.dlLastError(obj, val)
             
             obj.dlLastError = val;
             
        end
        
        function set.dlOutputLog(obj, val)
             
            obj.dlOutputLog = val;
             
        end
        
        function set.dlOutputs(obj, val)
             
            obj.dlOutputs = val;
             
        end
        
        function set.dlMexFuncName(obj, val)
             
            obj.dlMexFuncName = val;
             
        end
        
        function dlSave(obj)
        
            dlSaveFileNamePath = [obj.dlStudyDir, '/dlFile.mat'];
            p = load([obj.dlPath, '/params.mat']);
            save([obj.dlStudyDir, '/params.mat'], '-struct', 'p');
            save(dlSaveFileNamePath, 'obj');
            
        end
        
        function dlSaveCheckPoint(obj, dlCheckPointPath)
        
            p = load([obj.dlPath, '/params.mat']);
            save([obj.dlStudyDir, dlCheckPointPath, 'params.mat'], '-struct', 'p');
            fprintf("Checkpoint file saved in %s \n", [obj.dlStudyDir, dlCheckPointPath]);
            
        end
        
        function obj = dlLoadCheckPoint(obj, dlCheckPointPath)
            
            fprintf("Checkpoint file loaded from %s \n", [obj.dlStudyDir, dlCheckPointPath]);
            p = load([obj.dlStudyDir, dlCheckPointPath, 'params.mat']);
            save([obj.dlPath, '/params.mat'], '-struct', 'p');
            
        end
        
        function obj = dlLoad(obj, PathToFile)
            
            set(obj, 'dlPathToFile', [PathToFile, '/dlFile.mat']);
            o = load(obj.dlPathToFile);
            fprintf("DL object loaded from %s \n", PathToFile);
            obj = o.obj;
            
            fprintf("Params.mat file loaded from %s \n", PathToFile);
            p = load([PathToFile, '/solve/params.mat']);
            save([obj.dlPath, '/params.mat'], '-struct', 'p');
            obj.dlReInit();
            
        end
        
        function dlReInit(obj)
           
            [out, vars] = dsGetOutputList(obj.dlModel);
            set(obj, 'dlOutputs', out);
            set(obj, 'dlVariables', vars);
            obj.dlMexBridgeInit()
            
            obj.dlTrialNumber = 0;
            fprintf("\nReinitialized.\n");
            
        end
        
        function [s] = dlGetMexName(obj)
            
            obj.dlPath = [obj.dlStudyDir, '/solve'];
            addpath(obj.dlPath);
            d = dir(obj.dlPath);
            
            for i = 1:size(d, 1)
                if contains(d(i).name, 'mexw64')
                    s = d(i).name;
                    s = s(1:end-7);
                end
            end
            
        end
        
        function dlInit(obj, studydir) % Initializer TODO
            
            tspan = [0 10]; % Base time span for class construction and initialization.
            simulator_options = {'tspan', tspan, 'solver', 'rk1', 'dt', obj.dldT, ...
                        'downsample_factor', obj.dlDownSampleFactor, 'verbose_flag', 1, ...
                        'study_dir', studydir, 'mex_flag', 1};
            obj.dsData = dsSimulate(obj.dlModel, 'vary', [], simulator_options{:});
            
        end
        
        function dlMexBridgeInit(obj)
           
            set(obj, 'dlMexFuncName', obj.dlGetMexName());
            dsMexBridge('dlTempFunc.m', obj.dlMexFuncName);
            
        end
        
        function dlPlotAllPotentials(obj, mode)
           
            dlPotentialIndices = contains(obj.dlVariables, '_V');
            dlPotentialIndices(1) = 1;
            dlPotentials = obj.dlOutputs(dlPotentialIndices);
            dlLabels = obj.dlVariables(dlPotentialIndices);
            
            figure();
            patch([3 7 7 3], [-30 -30 +30 +30], [0.5 0.9 0.9]);hold("on");
            t = dlPotentials{1, 1};
            n = size(dlPotentials, 2);
            
            if strcmpi(mode, 'ifr')
                
                for i = 2:n

                    x = dlPotentials{1, i};
                    raster = computeRaster(t, x);
                    subplot(ceil(n/2), 2, i-1);

                    if size(raster, 1) > 0

                        pool = 1:size(x, 2);
                        O1 = 5e2 * NWepanechnikovKernelRegrRaster(t, raster, pool, 25, 1, 1);
                        plot(t, O1, 'o');

                    end

                    ylabel(dlLabels(i));

                end
                
                grid("on");title("iFR(s)");xlabel("time (ms)");
                fprintf("Done.\n");
                
            else
                
                for i = 2:n

                    x = dlPotentials{1, i};
                    subplot(ceil(n/2), 2, i-1);
                    plot(t, x);
                    ylabel(dlLabels(i));

                end
                
                grid("on");title("Voltage(s)");xlabel("time (ms)");
                fprintf("Done.\n");
                
            end

        end

        
        function dlSimulate(obj)
            
%             disp("Simulation ...");
            set(obj, 'dlOutputs', dlTempFunc(obj.dlOutputs));
%             disp("Done."); 
      
        end
        
        function kernel = dlTimeIntervalKernel(obj, dlTimeInterval)
            
            kernel = obj.dlOutputs{1, 1};
            kernel(kernel > dlTimeInterval(2)) = 0;
            kernel(kernel < dlTimeInterval(1)) = 0;
            kernel(kernel > 0) = 1; 
       
        end
        
        function out = dlApplyKernel(obj, dlOutput, dlKernel)
           
            xp = dlOutput;
            for j = 1:size(xp, 2)

                xp(:, j) = xp(:, j).*dlKernel;
                xp(xp == 0) = obj.dlBaseVoltage;

            end
            out = xp;
           
        end
        
        function out = dlApplyIFRKernel(obj, dlOutput)
              
            x = dlOutput;
            t = linspace(0, size(x, 1), size(x, 1))*obj.dldT*obj.dlDownSampleFactor;
            raster = computeRaster(t, x);

            if size(raster, 1) > 0

                pool = 1:size(x, 2);
                out = 5e2 * NWepanechnikovKernelRegrRaster(t, raster, pool, 25, 1, 1);
                
            else 

                out = zeros(1, size(x, 1));
                
            end
            
        end
        
        function out = dlApplyAverageFRKernel(obj, dlOutput)
           
            o1 = obj.dlApplyIFRKernel(dlOutput);
            out = mean(o1);
            
        end
        
        function dlCalculateOutputs(obj, dlOutputParameters)
           
            n = size(dlOutputParameters, 1);
            dlIndices = zeros(1, n);
            
            for i = 1:n
                
                dlIndices(i) = find(strcmpi(obj.dlVariables, dlOutputParameters{i, 1}));
                
            end
            
            set(obj, 'dlLastOutputs' , obj.dlOutputs(dlIndices));
            
            for i = 1:n
                
                dlTimeKernel = ceil(dlOutputParameters{i, 3}/(obj.dldT*obj.dlDownSampleFactor));
                dlOutputType = dlOutputParameters{i, 4};
                dlTempOutputs = obj.dlLastOutputs{i}(dlTimeKernel(1):dlTimeKernel(2), dlOutputParameters{i, 2});
         
                if strcmpi(dlOutputType, 'ifr')

                    obj.dlLastOutputs{i} = obj.dlApplyIFRKernel(dlTempOutputs);

                elseif strcmpi(dlOutputType, 'lfp')

                    obj.dlLastOutputs{i} = obj.dlApplyKernel(dlTempOutputs, dlTimeKernel);

                elseif strcmpi(dlOutputType, 'afr')

                    obj.dlLastOutputs{i} = obj.dlApplyAverageFRKernel(dlTempOutputs);

                else

                    fprintf("\tThis output type is not recognized. Trying to run ""%s.m"" ...\n", dlOutputType);
                    try

                        dsBridgeFunctionRun(dlOutputType);
                        obj.dlLastOutputs{i} = dlTempFuncBridge(dlTempOutputs);
                        fprintf("\t""%s.m"" runned succesfully for output no. %d\n", dlOutputType, i);

                    catch

                        fprintf("\tError occured. Function ""%s.m"" is not found or it contains error. check if you've created ""%s.m"" or entered correct output type.\n--->No valid output is calculated for this trial, response and error are going to be ""NaN""\n", dlOutputType, dlOutputType);
                        obj.dlLastOutputs{i} = "NaN";

                    end
                end
            end
        end
        
        function dlCalculateError(obj, dlTargetParams)
           
            n = size(dlTargetParams, 1);
            Error = 0;
            
            for i = 1:n
               
                TempError = 0;
                
                dlErrorType = dlTargetParams{i, 1};
                dlOutputIndices = dlTargetParams{i, 2};
                dlOutputTargets = dlTargetParams{i, 3};
                dlErrorWeight = dlTargetParams{i, 4};
                
                if strcmpi(dlErrorType, 'MAE')
                   
                    TempError = abs(obj.dlLastOutputs{dlOutputIndices} - dlOutputTargets);
                    
                elseif strcmpi(dlErrorType, 'MSE')
                    
                    TempError = abs(obj.dlLastOutputs{dlOutputIndices} - dlOutputTargets)^2;
                    
                elseif strcmpi(dlErrorType, 'Compare')
                    
                    x = dlOutputIndices;
                    c = 1;
                    
                    for j = dlOutputIndices
                        
                        x(c) = obj.dlLastOutputs{j};
                        
                        if c > 1
                            TempError = TempError + dlRampFunc(x(c) - x(c-1));
                        end
                        
                        c = c + 1;
                        
                    end
                
                elseif strcmpi(dlErrorType, 'Diff')
                    
                    j = dlOutputIndices;
                    TempError = abs(obj.dlLastOutputs{j(1)} - obj.dlLastOutputs{j(2)});
                    
                else
                    
                    fprintf("Undefined error type ""%s""\n", dlErrorType);
                    
                end
                
                Error = Error + TempError*dlErrorWeight;
                
            end
            
            obj.dlLastError = Error;
            obj.dlErrorsLog = [obj.dlErrorsLog, obj.dlLastError];
            
        end
        
        function out = dlAdaptiveLambda(obj)
            
            out = obj.dlLastLambda * obj.dlDeltaRatio;
            
        end
        
        function dlOutputGenerator(obj)
           
            n = size(obj.dlLastOutputs, 2);
            for i = 1:n
               
                obj.dlLastOutputs{i} = obj.dlLastOutputs{i} * rand(1) * 2;
                
            end
            
        end
        
        function dlRunSimulation(obj, dlVaryList, dlOutputParameters)
           
            fprintf("\n\tSingle trial running: \n");
            obj.dlUpdateParams(dlVaryList);
            obj.dlSimulate();
            obj.dlCalculateOutputs(dlOutputParameters);
            fprintf("\n\tSimulation outputs: ");
            disp(obj.dlLastOutputs);
            
        end
        
        function dlTrain(obj, dlInputParameters, dlOutputParameters, dlTargetParameters, dlTrainOptions) 
            
            try
               
                dlSimulationFlag = dlTrainOptions('dlSimulationFlag');
                
            catch
                
                dlSimulationFlag = 1;
                
            end
            
            try
               
                dlOfflineOutputGenerator = dlTrainOptions('dlOfflineOutputGenerator');
                
            catch
                
                dlOfflineOutputGenerator = 0;
                
            end
            
            try
               
                dlAdaptiveLambda = dlTrainOptions('dlAdaptiveLambda');
                
            catch
                
                dlAdaptiveLambda = 1;
                
            end
            
            try
               
                dlEpochs = dlTrainOptions('dlEpochs');
                
            catch
                
                dlEpochs = 10;
                fprintf("->Number of epochs was not determined in options map. Default dlEpochs = 10\n");
                
            end
            
            try
               
                dlBatchs = dlTrainOptions('dlBatchs');
                
            catch
                
                dlBatchs = size(dlInputParameters, 2);
                fprintf("->Batchs was not determined in options map, default dlBatchs = size(dlVaryList, 2)\n");
                
            end
            
            try
               
                dlLambda = dlTrainOptions('dlLambda');
                
            catch
                
                dlLambda = 0.001;
                fprintf("->Lambda was not determined in options map, default dlLambda = 1e-3\n");
                
            end
            
            try
               
                dlLearningRule = dlTrainOptions('dlLearningRule');
                
            catch
                
                dlLearningRule = 'DeltaRule';
                fprintf("->Learning rule was not determined in options map, default dlLearningRule = 'DeltaRule'\n");
                
            end
            
            try
               
                dlUpdateMode = dlTrainOptions('dlUpdateMode');
                
            catch
                
                dlUpdateMode = 'batch';
                fprintf("->Update mode was not determined in options map, default dlUpdateMode = 'batch'\n");
                
            end
            
            try
               
                dlOutputLogFlag = dlTrainOptions('dlOutputLogFlag');
                
            catch
                
                dlOutputLogFlag = 0;
                
            end
            
            try
               
                dlCheckpoint = dlTrainOptions('dlCheckpoint');
                
            catch
                
                dlCheckpoint = 'true';
                fprintf("->Checkpoint flag was not determined in options map, default dlCheckpoint = 'true'\n");
                
            end
            
            try
               
                dlCheckpointCoefficient = dlTrainOptions('dlCheckpointCoefficient');
                
            catch
                
                dlCheckpointCoefficient = 2.0;
                fprintf("->Checkpoint Coefficient for optimal state saving and loading was not determined in options map, default dlCheckpointCoefficient = 2\n");
                
            end
            
            if dlSimulationFlag ~= 1
               
                fprintf("->Simulation has been manually deactivated for this run.\n");
                if dlOfflineOutputGenerator == 1
               
                    fprintf("->Offline output generator (random) is activated. Outputs are only for debugging approaches. \n");
                
                end
                
            else
               
                if dlOfflineOutputGenerator == 1
               
                    fprintf("->Offline output generator (random) is activated but ignored as simulation is active. \n");
                
                end
                
            end
            
            if dlAdaptiveLambda == 1
               
                fprintf("->Adaptive lambda is active. Lambda (learning rate) will be changed based on volatility of model error.\n");
                
            end
            
            if dlOutputLogFlag == 1
               
                fprintf("->Outputs log will be saved.\n");
                
            end            
                
            for i = 1:dlEpochs
                
                fprintf("\tEpoch no. %d\n", i);
                for j = 1:dlBatchs
                
                    fprintf("\t\tBatch no. %d\t", j);
                    set(obj, 'dlTrialNumber', obj.dlTrialNumber + 1);
                    obj.dlUpdateParams(dlInputParameters{j});
                    
                    if dlSimulationFlag == 1
                        obj.dlSimulate();
                        obj.dlCalculateOutputs(dlOutputParameters);
                    elseif dlOfflineOutputGenerator == 1
                        obj.dlOutputGenerator();
                    end
                    
                    obj.dlCalculateError(dlTargetParameters{j});
                    fprintf("\tError = %f\n", obj.dlLastError);
                    
                    if dlOutputLogFlag
                        obj.dlOutputLog = [obj.dlOutputLog; obj.dlLastOutputs];
                    end
                    
                    if strcmpi(dlUpdateMode, 'trial')
                        
                        obj.dlUpdateError = obj.dlLastError;
                        if dlAdaptiveLambda == 1
                            dlLambda = obj.dlAdaptiveLambda();
                        end
                        obj.dlTrainStep(dlLearningRule, dlLambda);
                        
                    end
                    
                end
                
                dlAvgError = mean(obj.dlErrorsLog(end-2:end));
                fprintf("\t\tEpoch's Average Error = %f, Last lambda = %f\n", dlAvgError, dlLambda);
                
                if strcmpi(dlCheckpoint, 'true')
                    
                    if dlAvgError < obj.dlOptimalError

                        obj.dlOptimalError = dlAvgError;
                        obj.dlSaveOptimal();

                        if strcmpi(dlUpdateMode, 'batch')

                            obj.dlUpdateError = dlAvgError;
                            if dlAdaptiveLambda == 1
                                dlLambda = obj.dlAdaptiveLambda();
                            end
                            obj.dlTrainStep(dlLearningRule, dlLambda);

                        end

                    elseif dlAvgError > dlCheckpointCoefficient*obj.dlOptimalError

                        obj.dlLoadOptimal();

                    else

                        if strcmpi(dlUpdateMode, 'batch')

                            obj.dlUpdateError = dlAvgError;
                            if dlAdaptiveLambda == 1
                                dlLambda = obj.dlAdaptiveLambda();
                            end
                            obj.dlTrainStep(dlLearningRule, dlLambda);

                        end

                    end
                    
                else
                   
                    if strcmpi(dlUpdateMode, 'batch')

                        obj.dlUpdateError = dlAvgError;
                        if dlAdaptiveLambda == 1
                            dlLambda = obj.dlAdaptiveLambda();
                        end
                        obj.dlTrainStep(dlLearningRule, dlLambda);

                    end
                    
                end
            end
            
        end
        
        function dlTrainStep(obj, dlLearningRule, dlLambda)
           
            error = obj.dlUpdateError;
            obj.dlLastLambda = dlLambda;
            p = load([obj.dlPath, '/params.mat']);
            
            val = struct2cell(p.p);
            lab = fieldnames(p.p);
            l = find(contains(lab, '_netcon'));
            
            deltaL = 0;
            if strcmpi(dlLearningRule, 'DeltaRule')
            
                for i = l'

                    rng('shuffle');
                    w = val{i, 1};
                    delta = (randn(size(w)))*error*dlLambda;
                    wn = w + delta;
                    
                    wn(wn < 0) = 0;
                    wn(wn > 1) = 1;
                    val{i, 1} = wn;
                    deltaL = deltaL + sum(sum(abs(delta)));
                    
                end
                
            elseif strcmpi(dlLearningRule, 'BioDeltaRule')
            
                for i = l'

                    rng('shuffle');
                    w = val{i, 1};
                    delta = (1-w).*(randn(size(w)))*error*dlLambda;
                    wn = w + delta;
                    
                    wn(wn < 0) = 0;
                    wn(wn > 1) = 1;
                    val{i, 1} = wn;
                    deltaL = deltaL + sum(sum(abs(delta)));
                    
                end
                
            elseif strcmpi(dlLearningRule, 'RWDeltaRule')
            
                disp("TODO Rascorla-Wagner delta rule");
                
            else
                
                disp("TODO train step and learning 'else' part");
                disp(error);
                
            end
            
            if obj.dlLastDelta < 0
            
                obj.dlLastDelta = deltaL;
                    
            else
                
                obj.dlDeltaRatio = (obj.dlLastDelta / deltaL)^0.5;
                obj.dlLastDelta = deltaL;
                
            end
            
            q = cell2struct(val, lab);
            p.p = q;
            save([obj.dlPath, '/params.mat'], '-struct', 'p');
            
        end
        
        function c = dlGetConnectionsList(obj)
            
            p = load(obj.dlStudyDir + "/solve/params.mat");
            st = p.p;
            cl = fieldnames(st);
            c = [];
            
            for i = 1:size(cl, 1)
                if contains(cl(i), '_netcon')
                    c = [c; cl(i)];
                end
            end
            
        end
        
        function dlUpdateParams(obj, map) 
            
            fprintf("Updating parameters of %s", obj.dlPath);
            dsParamsModifier('dlTempFuncParamsChanger.m', map);
            dlTempFuncParamsChanger(obj.dlPath);
            fprintf("\tUpdated.\t"); 
            
        end
        
        function dlResetTraining(obj)
           
            obj.dlTrialNumber = 0;
            obj.dlOptimalError = 1e9;
            obj.dlOutputLog = [];
            obj.dlErrorsLog = [];
            
        end
        
        function dlSaveOptimal(obj)
           
            obj.dlSaveCheckPoint('/Optimal');
            
        end
        
        function dlLoadOptimal(obj)
            
            try
                obj.dlLoadCheckPoint('/Optimal');
            catch
                fprintf("--->No oprimal file exists. first run a training session with an active checkpoint flag.\n");
            end
            
        end
        
    end
end
