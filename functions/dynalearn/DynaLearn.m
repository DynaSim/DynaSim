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
        dlPathToFile = 'models/dlBaseModel';
        dlBaseVoltage = -77.4;
        dldT = .01; % Time step in ODEs (dt)
        
        dlDownSampleFactor = 10; 
        dlOptimalError = 1e9;
        dlUpdateError = 0;
        dlLastLambda = 1e-3;
        
        dlDeltaRatio = 1;
        dlLastDelta = -1;
        dlLambdaCap = 1e-2;
        dlSimulationTool = "mex";
        
    end
    
    methods

        function obj = DynaLearn(varargin) % Constructors, will be expanded
            
            fprintf("\n\n@DS.DL:Creating Dyna model object ... ");
            set(obj, 'dlPathToFile', 'models/dlBaseModel');
            
            if nargin == 0
                
                obj = obj.dlLoad(obj.dlPathToFile);
                
            elseif nargin == 1

                model_ = varargin{1};
                
                set(obj, 'dlModel', model_);
                obj.dlInit(obj.dlStudyDir);
                    
            elseif nargin == 2
                    
                model_ = varargin{1};  
                data_ = varargin{2};
                
                set(obj, 'dlModel', model_);
                set(obj, 'dlStudyDir', data_);
                
                obj.dlInit(obj.dlStudyDir);
                 
            elseif nargin == 3
                
                model_ = varargin{1};  
                data_ = varargin{2};
                mode_ = varargin{3};
                
                set(obj, 'dlModel', model_);
                set(obj, 'dlStudyDir', data_);
                set(obj, 'dlSimulationTool', mode_);
                
                if strcmpi(obj.dlSimulationTool, "mex")
                    obj.dlInit(obj.dlStudyDir);
                elseif strcmpi(obj.dlSimulationTool, "raw")
                    obj.dlRawInit(obj.dlStudyDir);
                else
                    fprintf("\n->Simulation tool ""%s"" is not valid. Try ""mex"" or ""raw"".\n ***No model created, try again.\n", mode_);
                    return
                end
                
            else
                fprintf("\n Invalid use of DynaNet; try one of the following ways:\n 1.Call with no inputs <> (for demo or reloading).\n 2.Call with <DynaSim struct>.\n 3.Call with <DynaSim struct, studydir>.\n 4.Call with <DynaSim Struct, studydir, simulation tool mode>.\n");
                return
            end
          
            [out, vars] = dsGetOutputList(obj.dlModel);
            set(obj, 'dlOutputs', out);
            set(obj, 'dlVariables', vars);
            
            if strcmpi(obj.dlSimulationTool, "mex")
                obj.dlMexBridgeInit();
            elseif strcmpi(obj.dlSimulationTool, "raw")
                obj.dlRawBridgeInit();
            end
            
            fprintf("\n@DS.DL:DynaLearn model created.\n");
            
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
        
        function set.dlSimulationTool(obj, val)
             
            obj.dlSimulationTool = val;
             
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
            
            set(obj, 'dlPathToFile', [PathToFile, '/dlFile.mat']);
            set(obj, 'dlPath', [PathToFile, '/solve']);
            
            fprintf("Params.mat file loaded from %s \n", PathToFile);
            p = load([PathToFile, '/solve/params.mat']);
            save([obj.dlPath, '/params.mat'], '-struct', 'p');
            obj.dlReInit();
            
        end
        
        function dlReInit(obj)
           
            [out, vars] = dsGetOutputList(obj.dlModel);
            set(obj, 'dlOutputs', out);
            set(obj, 'dlVariables', vars);
            
            if strcmpi("mex", obj.dlSimulationTool)
                obj.dlMexBridgeInit();
            elseif strcmpi("raw", obj.dlSimulationTool)
                obj.dlRawBridgeInit();
            end
            
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
        
        function [s] = dlGetRawName(obj)
            
            obj.dlPath = [obj.dlStudyDir, '/solve'];
            addpath(obj.dlPath);
            d = dir(obj.dlPath);
            
            for i = 1:size(d, 1)
                if contains(d(i).name, '.m')
                    s = d(i).name;
                    s = s(1:end-2);
                end
            end
            
        end
        
        function dlInit(obj, studydir) % Initializer with mex
            
            tspan = [0 10]; % Base time span for class construction and initialization.
            simulator_options = {'tspan', tspan, 'solver', 'rk1', 'dt', obj.dldT, ...
                        'downsample_factor', obj.dlDownSampleFactor, 'verbose_flag', 1, ...
                        'study_dir', studydir, 'mex_flag', 1};
            obj.dsData = dsSimulate(obj.dlModel, 'vary', [], simulator_options{:});
            
        end
        
        function dlRawInit(obj, studydir) % Initializer without mex
            
            tspan = [0 10]; % Base time span for class construction and initialization.
            simulator_options = {'tspan', tspan, 'solver', 'rk1', 'dt', obj.dldT, ...
                        'downsample_factor', obj.dlDownSampleFactor, 'verbose_flag', 1, ...
                        'study_dir', studydir, 'mex_flag', 0};
            obj.dsData = dsSimulate(obj.dlModel, 'vary', [], simulator_options{:});
            
        end
        
        function dlMexBridgeInit(obj)
           
            set(obj, 'dlMexFuncName', obj.dlGetMexName());
            dsMexBridge('dlTempFunc.m', obj.dlMexFuncName);
            
        end
        
        function dlRawBridgeInit(obj)
           
            set(obj, 'dlMexFuncName', obj.dlGetRawName());
            dsMexBridge('dlTempFunc.m', obj.dlMexFuncName);
            
        end
        
        function dlPlotAllPotentials(obj, mode, opts)
           
            dlPotentialIndices = contains(obj.dlVariables, '_V');
            dlPotentialIndices(1) = 1;
            dlPotentials = obj.dlOutputs(dlPotentialIndices);
            dlLabels = obj.dlVariables(dlPotentialIndices);
            

            t = dlPotentials{1, 1};
            n = size(dlPotentials, 2);
            m = ceil(n/6);
            
            for k = 1:m
                
                figure('Position', [0, 0, 1400, 700*(min(k*6, n-1) - (k-1)*6)/(6)]);
                
                if strcmpi(mode, 'ifr')

                    for i = (k-1)*6+1:min((k*6), n-1)

                        x = dlPotentials{1, i+1};
                        raster = computeRaster(t, x);
                        subplot((min(k*6, n-1) - (k-1)*6), 1, mod(i-1, (min(k*6, n-1) - (k-1)*6))+1);

                        if size(raster, 1) > 0

                            pool = 1:size(x, 2);
                            O1 = 5e2 * NWepanechnikovKernelRegrRaster(t, raster, pool, 49, 1, 1);
                            plot(t, O1, 'o');grid("on");

                        end

                        ylabel(dlLabels(i+1));

                    end

                    xlabel(mode + " in time (ms)");

                elseif strcmpi(mode, 'lfp')

                    for i = (k-1)*6+1:min((k*6), n-1)

                        x = dlPotentials{1, i+1};
                        subplot((min(k*6, n-1) - (k-1)*6), 1, mod(i-1, (min(k*6, n-1) - (k-1)*6))+1);
                        plot(t, x);grid("on");
                        ylabel(dlLabels(i+1));

                    end

                    xlabel(mode + " in time (ms)");
                    
                elseif strcmpi(mode, 'avglfp')

                    for i = (k-1)*6+1:min((k*6), 6)

                        x = dlPotentials{1, i+1};
                        subplot((min(k*6, n-1) - (k-1)*6), 1, mod(i-1, (min(k*6, n-1) - (k-1)*6))+1);
                        plot(t, mean(x, 2));grid("on");
                        ylabel(dlLabels(i+1));

                    end

                    disp("Temp edit for 6 subplots");
                    xlabel(mode + " in time (ms)");
                    return
                    
                elseif strcmpi(mode, 'avgfft')

                    dtf = ceil(1 / (obj.dldT*obj.dlDownSampleFactor));
                    
                    lf = opts("lf")*dtf;
                    hf = opts("hf")*dtf;
                    freqCap = 0;
                    
                    for i = (k-1)*6+1:min((k*6), 6)

                        x = dlPotentials{1, i+1};
                        fqs = linspace(1, 500, max(size(x)));
                        subplot((min(k*6, n-1) - (k-1)*6), 1, mod(i-1, (min(k*6, n-1) - (k-1)*6))+1);
                        ffts = abs(fft(mean(x, 2))) * min(size(x)) / 1000;
                        
                        plot(fqs(lf:hf), ffts(lf:hf));grid("on");
                        
                        if freqCap == 0
                            freqCap = max(ffts(lf:hf))*1.2;
                            ylim([0, freqCap]);
                        else
                            ylim([0, freqCap]);
                        end
                        
                        ylabel(dlLabels(i+1));

                    end

                    disp("Temp edit for 6 subplots; average fft");
                    xlabel(mode + " in frequency (Hz)");
                    return
                    
                else
                    
                    fprintf("--->Mode %s is not recognised. Try 'lfp' or other available options.\n", mode)
                    
                end
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
        
        function out = dlApplyIFRKernel(obj, dlOutput) % TODO: sampling rate change and adjustments
              
            x = dlOutput;
            t = linspace(0, size(x, 1), size(x, 1))*obj.dldT*obj.dlDownSampleFactor;
            raster = computeRaster(t, x);

            if size(raster, 1) > 0

                pool = 1:size(x, 2);
                out = 1.1e3 * NWepanechnikovKernelRegrRaster(t, raster, pool, 25, 1, 1);
                
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
            if ~ (out < obj.dlLambdaCap && out > 0)
                
                obj.dlDeltaRatio = 1;
                obj.dlLastDelta = -1;
                
                fprintf("--->Warning! Lambda (%f) exceeded the limit (%f). Consider divergence problems.\n", out, obj.dlLambdaCap);
                fprintf("--->To avoid problems, previous lambda will be used for the next update.\n");
                out = obj.dlLastLambda; 
                
            end
            
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
                try
                    obj.dlLambdaCap = dlTrainOptions('dlLambdaCap');
                catch
                    fprintf("-->Reminder: you have choosen adaptive lambda but forgot to determine a lambda cap. default dlLambdaCap = %f\n", obj.dlLambdaCap);
                end
            end
            
            if dlOutputLogFlag == 1
               
                fprintf("->Outputs log will be saved.\n");
                
            end            
                
            for i = 1:dlEpochs
                
                fprintf("\tEpoch no. %d (Total iterations for this model : %d)\n", i, obj.dlTrialNumber);
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
                    
                    if ~contains(lab{i, 1}, 'IO')
                        wn = w + delta;
                    end
                    
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
                    
                    if ~contains(lab{i, 1}, 'IO')
                        wn = w + delta;
                    end
                    
                    wn(wn < 0) = 0;
                    wn(wn > 1) = 1;
                    val{i, 1} = wn;
                    deltaL = deltaL + sum(sum(abs(delta)));
                    
                end

            %%% What follows are 4 Learning Rules (3 for E-cells and 1 for I-cells; WIP) from Clopath's paper:
            % Aljadeff et al. Cortical credit assignment by Hebbian, neuromodulatory and inhibitory plasticity. arXiv:1911.00307, 2019
            % In this study the final weight update is respectively:
            % wE(n+1) = min( [wE(n) + sum(delta_wE)]_+ , wE_max)
            % wI(n+1) = min( [wI(n) + sum(delta_wI)]_+ , wI_max)
            % Table 2 in the paper shows parameter values (PDF pg 21)
            % Model constraints can be seen paper's Table 1 (PDF pg 20)
            % Most important constraints (see further explanations below):
            % rho_NE < rho_ACh
            % 5A_ACh < A_NE
            % alpha_NE > alpha_ACh
            % alpha_ACh > alpha_Hebbian
            elseif strcmpi(dlLearningRule, 'ACh') % WIP
                % ACh learning rule for E-cells from Clopath's paper (LTP and LTD)
                % xE and Y sampled from Bernoulli distribution (they are 1 with probability f and 0 with probability 1-f)
                % xE: input (in Clopath's paper it was binary and stimulus specific)
                % yE: output (in Clopath's paper it was binary from the Heaviside step function)
                % f: reference spike prob for neuromodulation plasticity
                % eta: gating for ACh plasticity (input-output pairing, it was binary in Clopath's paper: 0/1)
                % alpha: learning rate for ACh plasticity
                % beta: LTD/LTP scaling factor
                % single expression for balanced LTD/LTP ('*' represents matrix multiplication):
                % delta = eta·alpha·(yE-f)*(xE - beta·<yE>/(1-<yE>)·(1-xE))

                % comments/things to consider:
                % - with f~0, no LTD even if y = 0
                % - the factor <yE>/(1-<yE>) ensures that for beta=1 potentiation and depression are balanced on average (not sure this is <yE> or f)
                % - pairing: if target is Y = 1, pairing prob is rho, sampling from that, eta will be 0/1, if Y = 0, then eta = 0
                % - M: if paired (eta=1) disinhibition (sampled rectified Gaussian distribution of mean and std = A)
            elseif strcmpi(dlLearningRule, 'NE') % WIP
                % NE learning rule for E-cells from Clopath's paper (only LTP and not stimulus specific)
                % xE sampled from Bernoulli distribution (they are 1 with probability f and 0 with probability 1-f)
                % xE: input (in Clopath's paper it was binary and stimulus specific)
                % yE: output (in Clopath's paper it was binary from the Heaviside step function)
                % f: reference spike prob for neuromodulation plasticity
                % eta: gating for NE plasticity (input-output pairing, it was binary in Clopath's paper: 0/1)
                % alpha: learning rate for NE plasticity
                % expression:
                % delta = eta·alpha·(yE-f)*xE

                % comments/things to consider:
                % - pairing: with rho prob indep of target Y
                % - M: if paired (eta=1) disinhibition (sampled rectified Gaussian distribution of mean and std = A)
            elseif strcmpi(dlLearningRule, 'Hebbian') % WIP
                % Hebbian learning rule for E-cells from Clopath's paper (not gated)
                % xE sampled from Bernoulli distribution (they are 1 with probability f and 0 with probability 1-f)
                % xE: input (in Clopath's paper it was binary and stimulus specific)
                % yE: output (in Clopath's paper it was binary from the Heaviside step function)
                % alpha: learning rate for Hebbian plasticity
                % expression:
                % delta = alpha·(yE-<yE>)*xE
            elseif strcmpi(dlLearningRule, 'Inhibitory') % WIP
                % learning rule for I-cells from Clopath's paper (based on detailed E-I balance)
                % xE and xI sampled from Bernoulli distribution (they are 1 with probability f and 0 with probability 1-f)
                % xE: input (in Clopath's paper it was binary and stimulus specific)
                % xI: input (in Clopath's paper it was binary and stimulus specific)
                % excitatory current iE = wE*xE
                % inhibitory current iI = wI*xI
                % E/I balance line (reference): iI = a·iE + b (with a < 1)
                % alpha: learning rate for Inhibitory plasticity
                % expression:
                % delta = alpha·((a·iE + b)-iI)*xI

            elseif strcmpi(dlLearningRule, 'UncertaintyReduction')

                dlLambdaCap = 1; % allowing learning rates to be in [0,1]
                dlAdaptiveLambda = 0; % disabling adaptive lambda as this learning rule controls the lambda values through uncertainty reduction
                uncBaseline = 0.5; % reference value for uncertainty reduction (point of maximum uncertainty)
                scalingFactor = max([uncBaseline, 1-uncBaseline]); % used to keep uncReduct in [0,1]
                stochasticFactor = 0.2; % stochastic modulation

                dlLambda0 = dlLambda;
                for i = l'
                    rng('shuffle'); %% TODO we shouldn't shuffle all the time
                    w = val{i, 1};

                    % lambda update from previous w
                    if isscalar(dlLambda0) % first time
                        mu = dlLambda0; %% TODO add mu as a dl parameter?
                        lambda = dlLambda0*ones(size(w));
                    else % subsequent times
                        lambda = dlLambda{i, 1};
                    end
                    uncReduct = abs(w-uncBaseline)/scalingFactor; % Uncertainty reduction in [0,1]
                    % adapting lambda based on uncertainty reduction (the lower the uncertainty, the faster it adapts)
                    lambda = lambda + mu*(uncReduct-lambda); % adapting lambda based on uncertainty reduction

                    % w update based on the new lambda
                    delta = (1 + stochasticFactor*randn(size(w))).*lambda.*(1-w)*error; % stochastic delta
                    wn = w + delta;

                    % rectifying values that are out of the [0,1] bounds
                    wn(wn < 0) = 0;
                    wn(wn > 1) = 1;

                    % saving
                    val{i, 1} = wn;
                    dlLambda{i, 1} = lambda;

                    deltaL = deltaL + sum(sum(abs(delta)));
                end

            elseif strcmpi(dlLearningRule, 'RWDeltaRule')
            
                disp("TODO Rascorla-Wagner delta rule");
                
            elseif strcmpi(dlLearningRule, 'NewRule')
            
                disp("TODO new rule");
                
            else
                
                disp("TODO train step and learning 'else' part");
                
            end
            
            if obj.dlLastDelta < 0
            
                obj.dlLastDelta = deltaL;
                    
            else
                
                obj.dlDeltaRatio = min((obj.dlLastDelta / deltaL)^0.5, 2);
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
        
        function dlPlotErrors(obj)
           
            figure('position', [0, 0, 1400, 700]);
            plot(obj.dlErrorsLog);grid("on");
            title("Errors in trials");
            
        end
        
        function dlPlotBatchErrors(obj, dlBatchs)
           
            figure('position', [0, 0, 1400, 700]);
            n = max(size(obj.dlErrorsLog));
            x = zeros(1, ceil(n/dlBatchs));
            
            for i = 0:dlBatchs:n-dlBatchs
                x(ceil((i+1)/dlBatchs)) = mean(obj.dlErrorsLog(i+1:i+dlBatchs));
            end
            
            plot(x);grid("on");
            title("Errors in batchs");
            
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
                fprintf("--->No oprimal file exists. first run a training session with an active checkpoint flag to save an optimal checkpoint.\n");
            end
            
        end
        
    end
end
