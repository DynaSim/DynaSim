classdef DynaLearn < matlab.mixin.SetGet

    %

    % DynaLearn on DynaSim; Stable

    % Integrated class for interactive and trainable DynaSim network models

    % Toolbox is under developement for the future versions.

    %
    
    properties
        
        dlModel = []; % Dynasim's model struct, containing model's populations (layers) and connections, mechanisms ... .
        dlStudyDir = 'dlModels/'; % Output data of last simulation.
        dlErrorsLog = []; % Log of errors since first trial.
        dlConnections = []; % List of connection names.
    
        dlTrialNumber = 0; % Last trial; how many times this model have been trained with stimuli.
        dlLastOutputs = []; % Last output (eg. spike vector) of this model.
        dlLastError = 0; % Last error of this model for target, equivalent to the last value of errors_log property
        dlOutputLog = {}; % Last expected output that this model should've generated.
        
        dsData = []; % Last simulation outputs
        dlOutputs = []; % Mex outputs
        dlVariables = []; % Mex variable labels
        dlMexFuncName = []; % Name of Mex function (e.g **********_mex.mex64

        dlCurrentSessionValidTrials = 0;
        dlEnhancedDeltaRuleState = 1;
        dlLastWeightChanges = []; % To keep track of optimal changes
        dlExcludeDiverge = 0;

        dlLastErrorsLog = 1;
        dlLastCustomLog = 1;
        dlWeightsValues = []; % Weights values history {[Npre,Npost,1+Epochs]}
        dlWeightsVariables = []; % Weights variables
        
        dlPath = []; % Path which contains params.mat, mexfuncs, solver ...
        dlPathToFile = 'models/dlBaseModel';
        dlBaseVoltage = -77.4;
        dldT = .01; % Time step in ODEs (dt)
        
        dlDownSampleFactor = 10; % dS parameter for downsampling computations
        dlOptimalError = 1e9; % dL optimal error of training, initially it is just an irrelevant high number
        dlLastOptimalTrial = 1; % The trial witl best results
        dlUpdateError = 0; % The error which is used to update last state
        
        dlLastLambda = 1e-7; % Last lambda parameter
        dlDeltaRatio = 1;
        dlLastDelta = -1;
        dlLambdaCap = 1;

        dlLastOutputLog = {};
        dlOptimalWeightChanges = [];
        dlOptimalWeightChangesFlag = 0;
        dlModelName = "dlModelName";

        % MetaLR
        dlMetaMu = 0.01;
        dlMetaLR = 0.01;
        error0 = nan;
        
        dlGraph = []; 
        dlCustomLog = {};
        dlCustomLogLabel = [];
        dlSimulationTool = "mex";

        dlMaxFrequency = 500; % Can be changed
        dlParams = [];

    end

    methods (Static)

        function obj = dlLoader(PathToFile)
            
            o = load([PathToFile, '/dlFile.mat']);
            obj = o.obj;
            set(obj, 'dlPathToFile', [PathToFile, '/dlFile.mat']);
            fprintf("DL object loaded from %s \n", PathToFile);
            
            set(obj, 'dlPathToFile', [PathToFile, '/dlFile.mat']);
            set(obj, 'dlPath', [PathToFile, '/solve']);
            
            fprintf("Params.mat file loaded from %s \n", PathToFile);
            p = load([PathToFile, '/params.mat']);
            save([obj.dlPath, '/params.mat'], '-struct', 'p');
            obj.dlReInit();
            
        end

    end

    methods

        function obj = DynaLearn(varargin) % Constructors, will be expanded
            
            fprintf("\n\n@DS.DL:Creating DynaLearn model object ... \n");
            set(obj, 'dlPathToFile', 'models/dlBaseModel');
            
            if nargin == 0
                
                return
                
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
                set(obj, 'dlPath', data_);

                if strcmpi(obj.dlSimulationTool, "mex")
                    obj.dlInit(obj.dlStudyDir);
                elseif strcmpi(obj.dlSimulationTool, "raw")
                    obj.dlRawInit(obj.dlStudyDir);
                else
                    fprintf("\n->Simulation tool ""%s"" is not valid. Try ""mex"" or ""raw"".\n ***No model created, try again.\n", mode_);
                    return
                end
            
            elseif nargin == 4
                
                model_ = varargin{1};  
                data_ = varargin{2};
                mode_ = varargin{3};
                name_ = varargin{4};
                
                set(obj, 'dlModel', model_);
                set(obj, 'dlStudyDir', data_);
                set(obj, 'dlSimulationTool', mode_);
                set(obj, 'dlPath', data_);

                set(obj, 'dlModelName', name_);

                if strcmpi(obj.dlSimulationTool, "mex")
                    obj.dlInit(obj.dlStudyDir);
                elseif strcmpi(obj.dlSimulationTool, "raw")
                    obj.dlRawInit(obj.dlStudyDir);
                else
                    fprintf("\n->Simulation tool ""%s"" is not valid. Try ""mex"" or ""raw"".\n ***No model created, try again.\n", mode_);
                    return
                end

            else

                fprintf("\n Invalid use of DynaLearn; try one of the following ways:\n 1.Call with no inputs <> (for demo or reloading).\n 2.Call with <DynaSim struct>.\n 3.Call with <DynaSim struct, studydir>.\n 4.Call with <DynaSim Struct, studydir, simulation tool mode>.\n");
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
            
            obj.dlGraphConstructor();
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

        function set.dlModelName(obj, val)
            
             if ~strcmpi(class(val), 'string') && ~ischar(val)
                error('Model name a string');
             end
             obj.dlModelName = val;
             
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

            fprintf("\n->Model saved in ""%s"".\n", dlSaveFileNamePath);
            
        end
        
        function dlSaveCheckPoint(obj, dlCheckPointPath)
        
            fprintf("\t\tCheckpoint file saved in %s \n", [obj.dlStudyDir, dlCheckPointPath]);
            p = load([obj.dlPath, '/params.mat']);
            save([obj.dlStudyDir, dlCheckPointPath, 'params.mat'], '-struct', 'p');
            save([obj.dlStudyDir, dlCheckPointPath, 'object.mat'], 'obj');
            
        end
        
        function out = dlLoadCheckPoint(obj, dlCheckPointPath)
            
            fprintf("\t\tCheckpoint file loaded from %s \n", [obj.dlStudyDir, dlCheckPointPath]);

            try
                dlObj = load([obj.dlStudyDir, dlCheckPointPath, 'object.mat']);
                out = dlObj.obj;
            catch
                disp("---->Optimal skipped (not loaded)");
                out = obj;
            end
            
            if obj.dlExcludeDiverge == 1

                obj.dlErrorsLog = obj.dlLastErrorsLog;
                obj.dlCustomLog = obj.dlLastCustomLog;
                save([obj.dlStudyDir, dlCheckPointPath, 'object.mat'], 'obj');

            else

                obj.dlErrorsLog = obj.dlErrorsLog;
                obj.dlCustomLog = obj.dlCustomLog;

            end

            disp([obj.dlStudyDir, dlCheckPointPath, 'params.mat']);

            try

                p = load([obj.dlStudyDir, dlCheckPointPath, 'params.mat']);
                obj.dlParams = p.p;
                save([obj.dlPath, '/params.mat'], '-struct', 'p');

            catch

                disp("---->Optimal skipped (not loaded)");  
                p = load([obj.dlPath, '/params.mat']);
                obj.dlParams = p.p;
                % save([obj.dlPath, '/params.mat'], '-struct', 'p');

            end
            
        end

        function out = dlLoadCheckPointLX(obj, dlCheckPointPath)
            
            fprintf("\t\tCheckpoint file loaded from %s \n", [obj.dlStudyDir, dlCheckPointPath]);
            dlObj = load([obj.dlStudyDir, dlCheckPointPath, 'object.mat']);
            out = dlObj.obj;
            
            if obj.dlExcludeDiverge == 1

                k = size(obj.dlLastCustomLog, 2);

                obj.dlLastErrorsLog = [obj.dlLastErrorsLog, obj.dlOptimalError*1.001];
               
                try
                    obj.dlLastCustomLog(:, k+1) = obj.dlLastCustomLog(:, k);
                catch

                end

                obj.dlErrorsLog = obj.dlLastErrorsLog;
                obj.dlCustomLog = obj.dlLastCustomLog;
                
                save([obj.dlStudyDir, dlCheckPointPath, 'object.mat'], 'obj');

            else

                obj.dlErrorsLog = obj.dlErrorsLog;
                obj.dlCustomLog = obj.dlCustomLog;

            end

            p = load([obj.dlStudyDir, dlCheckPointPath, 'params.mat']);
            obj.dlParams = p.p;
            save([obj.dlPath, '/params.mat'], '-struct', 'p');
            
        end
        
        function y = dlLoad(obj, PathToFile)
            
            set(obj, 'dlPathToFile', [PathToFile, '/dlFile.mat']);
            o = load(obj.dlPathToFile);
            fprintf("DL object loaded from %s \n", PathToFile);
            y = o.obj;
            
            set(obj, 'dlPathToFile', [PathToFile, '/dlFile.mat']);
            set(obj, 'dlPath', [PathToFile, '/solve']);
            
            fprintf("Params.mat file loaded from %s \n", PathToFile);
            p = load([PathToFile, '/solve/params.mat']);
            save([obj.dlPath, '/params.mat'], '-struct', 'p');
            obj.dlReInit();

            obj.dlParams = p.p;
            
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
            obj.dlGraphConstructor();

            try

                p = load([obj.dlPath, '/params.mat']);
                obj.dlParams = p.p;

            catch

                fprintf("->Warning: No params.mat file found. Run a simulation via obj.dlSimulate() or train the model for few trials (e.g 1) to establish one.\n");
            end

            fprintf("\nReinitialized.\n");
            
        end
        
        function [s] = dlGetMexName(obj)
            
            obj.dlPath = [obj.dlStudyDir, '/solve'];
            addpath(obj.dlPath);
            d = dir(obj.dlPath);
            s = "";
            
            for i = 1:size(d, 1)

                if contains(d(i).name, 'mexmaci64')

                    s = d(i).name;
                    s = s(1:end-10);

                elseif contains(d(i).name, 'mexw64')

                    s = d(i).name;
                    s = s(1:end-7);

                elseif contains(d(i).name, 'mexa64')

                    s = d(i).name;
                    s = s(1:end-7);

                end

            end
            
        end
        
        function [s] = dlGetRawName(obj)
            
            obj.dlPath = [obj.dlStudyDir, '/solve'];
            addpath(obj.dlPath);
            d = dir(obj.dlPath);
            s = "";
            
            for i = 1:size(d, 1)
                if contains(d(i).name, '.m')
                    s = d(i).name;
                    s = s(1:end-2);
                end
            end
            
        end
        
        function dlInit(obj, studydir) % Initializer with mex
            
            tspan = [0 1000]; % Base time span for class construction and initialization.
            simulator_options = {'tspan', tspan, 'solver', 'rk1', 'dt', obj.dldT, ...
                        'downsample_factor', obj.dlDownSampleFactor, 'verbose_flag', 1, ...
                        'study_dir', studydir, 'mex_flag', 1, 'mex_dir', obj.dlPath};
            obj.dsData = dsSimulate(obj.dlModel, 'vary', [], simulator_options{:});
            
        end
        
        function dlRawInit(obj, studydir) % Initializer without mex
            
            tspan = [0 1000]; % Base time span for class construction and initialization.
            simulator_options = {'tspan', tspan, 'solver', 'rk1', 'dt', obj.dldT, ...
                        'downsample_factor', obj.dlDownSampleFactor, 'verbose_flag', 1, ...
                        'study_dir', studydir, 'mex_flag', 0, 'mex_dir', obj.dlPath};
            obj.dsData = dsSimulate(obj.dlModel, 'vary', [], simulator_options{:});
            
            try

                p = load([obj.dlPath, '/params.mat']);
                obj.dlParams = p.p;

            catch

                fprintf("->Warning: No params.mat file found. Run a simulation via obj.dlSimulate() or train the model for few trials (e.g 1) to establish one.\n");
            
            end

        end
        
        function dlMexBridgeInit(obj)
           
            set(obj, 'dlMexFuncName', obj.dlGetMexName());
            dsMexBridge('dlTempFunc.m', obj.dlMexFuncName);
            
        end
        
        function dlRawBridgeInit(obj)
           
            set(obj, 'dlMexFuncName', obj.dlGetRawName());
            dsMexBridge('dlTempFunc.m', obj.dlMexFuncName);
            
        end

        function y = dlPopulationLabelTrim(~, label)

            if iscell(label)
                label = label{1};
            end

            indx = strfind(label, "_");
            N = length(indx);
            indx = indx(N - 1);
            y = label(1:indx-1);

        end
        
        function dlPlotAllPotentials(obj, mode, opts)
           
            dlPotentialIndices = endsWith(obj.dlVariables, '_V');
            dlPotentialIndices(1) = 1;
            dlPotentials = obj.dlOutputs(dlPotentialIndices);
            dlLabels = obj.dlVariables(dlPotentialIndices);

            t = dlPotentials{1, 1};
            n = size(dlPotentials, 2);
            q = 16;

            m = ceil(n/q);
            title_ = "Model: " + obj.dlModelName;
            
            for k = 1:m
                
                if strcmpi(mode, 'ifr')

                    figure('Position', [0, 0, 1700, 700*(min(k*7, n-1) - (k-1)*7)/(7)]);
                
                    for i = (k-1)*q+1:min((k*q), n-1)

                        x = dlPotentials{1, i+1};
                        raster = computeRaster(t, x);
                        subplot((min(k*q, n-1) - (k-1)*q), 1, mod(i-1, (min(k*q, n-1) - (k-1)*q))+1);

                        if size(raster, 1) > 0

                            pool = 1:size(x, 2);
                            O1 = 5e2 * dlNWRasterToIFR(t, raster, pool, 49, 1, 1);
                            plot(t, O1, 'o');grid("on");

                        end

                        ylabel(obj.dlPopulationLabelTrim(dlLabels(i+1)));
                        title(title_);

                    end

                    xlabel(mode + " in time (ms)");

                elseif strcmpi(mode, 'lfp')

                    figure('Position', [0, 0, 1700, 700*(min(k*q, n-1) - (k-1)*q)/(q)]);
                
                    for i = (k-1)*q+1:min((k*q), n-1)

                        x = dlPotentials{1, i+1};
                        subplot((min(k*q, n-1) - (k-1)*q), 1, mod(i-1, (min(k*q, n-1) - (k-1)*q))+1);
                        plot(t, x);
                        grid("on");

                        try

                            ylabel(obj.dlPopulationLabelTrim(dlLabels(i+1)));

                        catch

                            ylabel(dlLabels(i+1));

                        end

                        title(title_);

                    end

                    xlabel(mode + " in time (ms)");
                    
                elseif strcmpi(mode, 'avglfp')

                    figure('Position', [0, 0, 1700, 700*(min(k*7, n-1) - (k-1)*7)/(7)]);
                
                    for i = (k-1)*q+1:min((k*q), n-1)

                        x = dlPotentials{1, i+1};
                        subplot((min(k*q, n-1) - (k-1)*q), 1, mod(i-1, (min(k*q, n-1) - (k-1)*q))+1);
                        plot(t, mean(x, 2));
                        grid("on");

                        ylabel(obj.dlPopulationLabelTrim(dlLabels(i+1)));
                        title(title_);

                    end

                    disp("Temp edit for q subplots");
                    xlabel(mode + " in time (ms)");
                    return
                    
                elseif strcmpi(mode, 'avgfft')

                    figure('Position', [0, 0, 1700, 700*(min(k*7, n-1) - (k-1)*7)/(7)]);
                
                    dtf = ceil(1 / (obj.dldT*obj.dlDownSampleFactor));
                    
                    lf = opts("lf")*dtf;
                    hf = opts("hf")*dtf; 
                    freqCap = 0;
                    
                    for i = (k-1)*q+1:min((k*q)-1, n-1)

                        x = dlPotentials{1, i+1};
                        fqs = linspace(1, 500, max(size(x)));
                        subplot((min(k*q, n-1) - (k-1)*q), 1, mod(i-1, (min(k*q, n-1) - (k-1)*q))+1);
                        ffts = abs(fft(mean(x, 2))) * min(size(x)) / 1000;

                        try

                            yf = smooth(ffts(lf:hf));
                            area(fqs(lf:hf), yf);
                            grid("on");

                            if freqCap == 0
                                freqCap = max(ffts(lf:hf))*1.2;
                                ylim([0, freqCap]);
                            else
                                ylim([0, freqCap]);
                            end

                            ylabel(obj.dlPopulationLabelTrim(dlLabels(i+1)));
                            title(title_);

                        catch

                           fprintf(" Error in plotter* "); 

                        end

                    end

                    disp("Temp edit for q subplots; average fft");
                    xlabel(mode + " in frequency (Hz)");
                    return
                 
                elseif strcmpi(mode, 'raster')

                    imgRaster = zeros(1, length(t));
                    yticklabels_ = string(length(dlLabels)-1);
                    yticks_ = zeros(1, length(dlLabels)-1);
                    yticks_l = zeros(1, length(dlLabels)-1);

                    for i = (k-1)*q+1:min((k*q)-1, n-1)

                        x = dlPotentials{1, i+1};

                        yticks_l(i) = size(imgRaster, 1);

                        if size(imgRaster, 1) == 1

                            imgRaster(1:end+size(x, 2)-1, :) = x';

                        else

                            imgRaster(end+1:end+size(x, 2), :) = x';

                        end

                        yticklabels_(i) = obj.dlPopulationLabelTrim(dlLabels(i+1));
                        yticks_(i) = size(imgRaster, 1);

                    end

                    figure('Position', [0, 0, 1700, 1400]);
                    imagesc(1-imgRaster>4, "XData", t);
                    colormap("gray");
                    xn = size(imgRaster, 2);

                    hold on;

                    for i = (k-1)*q+1:min((k*q), n-1)

                        
                        line([0, xn], [yticks_(i)+.5, yticks_(i)+.5], 'Color', 'r');

                    end

                    yticks((yticks_+yticks_l)/2);
                    yticklabels(yticklabels_);
                    colormap("gray");

                    ylim([.5 yticks_(min((k*q), n-1))+.5]);
                    xlabel(mode + " in time (ms)");
                    title(title_);

                elseif strcmpi(mode, 'lfpmap')

                    imgRaster = zeros(1, length(t));
                    yticklabels_ = string(length(dlLabels)-1);
                    yticks_ = zeros(1, length(dlLabels)-1);
                    yticks_l = zeros(1, length(dlLabels)-1);

                    for i = (k-1)*q+1:min((k*q)-1, n-1)

                        x = dlPotentials{1, i+1};

                        yticks_l(i) = size(imgRaster, 1)*100;

                        if size(imgRaster, 1) == 1

                            imgRaster(1:end+size(x, 2)-1, :) = x';

                        else

                            imgRaster(end+1:end+size(x, 2), :) = x';

                        end

                        yticklabels_(i) = obj.dlPopulationLabelTrim(dlLabels(i+1));
                        yticks_(i) = (size(imgRaster, 1)+1)*100;

                    end

                    figure('Position', [0, 0, 1700, 1400]);
                    xn = max(t);

                    for i = 1:size(imgRaster, 1)

                        plot(t, -imgRaster(i, :)+i*100);
                        hold("on");

                    end
                                    
                    for i = (k-1)*q+1:min((k*q), n-1)

                        
                        line([0, xn], [yticks_(i), yticks_(i)], 'Color', 'black','LineStyle','--');

                    end

                    % ylim([.5 yticks_(min((k*q), n-1))+.5]);
                    disp(yticks_);
                    yticks((yticks_+yticks_l)/2);
                    yticklabels(yticklabels_);
                    xlabel(mode + " in time (ms)");
                    title(title_);
                    set(gca, 'YDir','reverse');

                elseif strcmpi(mode, 'lfpsave')

                    imgRaster = zeros(1, length(t));
                    for i = (k-1)*q+1:min((k*q), n-1)

                        x = dlPotentials{1, i+1};
                        imgRaster(end+1:end+size(x, 2), :) = x';
                        imgRaster(end+1, :) = 5;

                    end

                    save(obj.dlPath + "/signals_" + opts.name, "imgRaster");

                else
                    
                    fprintf("--->Mode %s is not recognised. Try 'lfp' or other available options.\n", mode)
                    
                end
            end
        end

        
        function dlSimulate(obj)
            
            obj.dlOutputs = cell(1,numel(obj.dlVariables));
            set(obj, 'dlOutputs', dlTempFunc(obj.dlOutputs));
            p = load([obj.dlStudyDir, '/solve/params.mat']);
            obj.dlParams = p.p;

            disp("Done."); 
      
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
        
        function out = dlMeanFR(obj, dlOutput)
            
            x = dlOutput;
            t = linspace(0, size(x, 1), size(x, 1))*obj.dldT*obj.dlDownSampleFactor;
            raster = computeRaster(t, x);

            if size(raster, 1) > 0
                nPop = size(x,2);
                out = size(raster,1)/(1e-3*(t(end)-t(1)))/nPop;
            else
                out = 0;
            end

        end

        function out = dlApplyIFRKernel(obj, dlOutput)
              
            x = dlOutput;

            for i = 1:size(x, 2)

                x(:, i) = x(:, i) - min(x(:, i));
                x(:, i) = x(:, i) / 40;
                x(:, i) = x(:, i)*90 - 70;

            end

            t = linspace(0, size(x, 1), size(x, 1))*obj.dldT*obj.dlDownSampleFactor;
            tScale = .9;% / (size(x, 1)*obj.dldT*obj.dlDownSampleFactor/1000);
            raster = computeRaster(t, x);
            pool = 1:size(x, 2);

            if size(raster, 1) > 0

                out = tScale * 1e3 * dlNWRasterToIFR(t, raster, pool, 25, 1, 1);
                
            else 

                out = zeros(1, size(x, 1));
                
            end

            % disp(tScale*out);
            
        end
        
        function out = dlApplyAverageFRKernel(obj, dlOutput)

            % o1 = obj.dlApplyIFRKernel(dlOutput);
            o1 = obj.dlMeanFR(dlOutput);
            out = mean(o1);
            
        end

        function out = dlApplyNaturalFrequencyKernel(obj, dlOutput)

            out = dlCalcNaturalFrequency(dlOutput, obj.dldT, obj.dlDownSampleFactor);

        end
        
        function dlCalculateOutputs(obj, dlOutputParameters)
           
            n = size(dlOutputParameters, 1);
            dlIndices = zeros(1, n);
            
            % disp(obj.dlVariables);
            % for i = 1:n
            %     disp(dlOutputParameters{i, 1});
            % end

            for i = 1:n
                
                try

                    dlIndices(i) = find(strcmpi(obj.dlVariables, dlOutputParameters{i, 1}));

                catch

                    fprintf("\n-->Check your model or output parameters, there is a problem about their name. Session is going to be invalid.\n");
                    error("\n@ds.dl: Model parameters do not match outputs or its variables.\n");

                end
                
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
                    
                elseif strcmpi(dlOutputType, 'astd')

                    obj.dlLastOutputs{i} = mean(std(dlTempOutputs), 'all');

                elseif strcmpi(dlOutputType, 'av')

                    obj.dlLastOutputs{i} = mean(mean(dlTempOutputs), 'all');

                elseif strcmpi(dlOutputType, 'afr')

                    obj.dlLastOutputs{i} = obj.dlApplyAverageFRKernel(dlTempOutputs);

                elseif strcmpi(dlOutputType, 'fnat')

                    obj.dlLastOutputs{i} = obj.dlApplyNaturalFrequencyKernel(dlTempOutputs);

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

                try

                    dlLowerFreq = dlTargetParams{i, 5};
                    dlUpperFreq = dlTargetParams{i, 6};

                catch

                    dlLowerFreq = 1;
                    dlUpperFreq = 100;

                end

                try

                    dlLowerFreq1 = dlTargetParams{i, 5};
                    dlUpperFreq1 = dlTargetParams{i, 6};
                    dlLowerFreq2 = dlTargetParams{i, 7};
                    dlUpperFreq2 = dlTargetParams{i, 8};

                catch

                    dlLowerFreq1 = 1;
                    dlUpperFreq1 = 49;
                    dlLowerFreq2 = 50;
                    dlUpperFreq2 = 100;

                end
                
                if strcmpi(dlErrorType, 'MAE')
                   
                    TempError = abs(obj.dlLastOutputs{dlOutputIndices} - dlOutputTargets);
                    
                elseif strcmpi(dlErrorType, 'Error') % Raw error keeping sign

                    TempError = obj.dlLastOutputs{dlOutputIndices} - dlOutputTargets;

                elseif strcmpi(dlErrorType, 'MSE')
                    
                    TempError = abs(obj.dlLastOutputs{dlOutputIndices} - dlOutputTargets)^2;
                
                elseif strcmpi(dlErrorType, 'MQE')
                    
                    TempError = abs(obj.dlLastOutputs{dlOutputIndices} - dlOutputTargets)^4;
                
                elseif strcmpi(dlErrorType, 'Compare')
                    
                    x = dlOutputIndices;
                    c = 1;
                    
                    for j = dlOutputIndices
                        
                        x(c) = obj.dlLastOutputs{j};

                        if c > 1
                            TempError = TempError + dlRampFunc(x(c) - x(c-1)).^2;
                        end
                        
                        c = c + 1;
                        
                    end
                
                elseif strcmpi(dlErrorType, 'Diff')
                    
                    j = dlOutputIndices;
                    TempError = abs(obj.dlLastOutputs{j(1)} - obj.dlLastOutputs{j(2)});
                
                elseif strcmpi(dlErrorType, 'TotalSpikesPenalty')
                    
                    TempError = 0;

                    for j = dlOutputIndices

                      TempError = squeeze(TempError  + obj.dlLastOutputs{j});
        
                    end

                    fprintf(" rEp=%d ", TempError);
                    TempError = abs(TempError - dlOutputTargets)^2;

                elseif strcmpi(dlErrorType, 'EPenalty')
                    
                    argsPSR = struct();

                    argsPSR.lf1 = dlLowerFreq;
                    argsPSR.hf1 = dlUpperFreq;

                    TempError = mean(dlEPowerSpectrum(obj, argsPSR), 'all');
                    fprintf(" gEp=%d ", TempError);
                    TempError = abs(TempError - dlOutputTargets)^2;
                    fprintf(" dlTarg=%d ", dlOutputTargets);
               
                elseif strcmpi(dlErrorType, 'RPenalty')
                    
                    argsPSR = struct();

                    argsPSR.lf1 = dlLowerFreq1; 
                    argsPSR.hf1 = dlUpperFreq1;
                    argsPSR.lf2 = dlLowerFreq2;
                    argsPSR.hf2 = dlUpperFreq2;

                    TempError = mean(dlRPowerSpectrum(obj, argsPSR), 'all');
                    fprintf(" gRp=%d f:[%d-%d Hz / %d-%d Hz] ", TempError, dlLowerFreq1, dlUpperFreq1, dlLowerFreq2, dlUpperFreq2);
                    TempError = abs(TempError - dlOutputTargets)^2;
                    fprintf(" dlTarg=%d ", dlOutputTargets);

                else
                    
                    fprintf("Undefined error type ""%s""\n", dlErrorType);
                    
                end
                
                Error = Error + TempError*dlErrorWeight;
                
            end
            
            if isnan(Error)

                Error = 1e6;

            end

            obj.dlLastError = Error;
            obj.dlErrorsLog = [obj.dlErrorsLog, obj.dlLastError];
            
        end
        
        function out = dlAdaptiveLambda(obj)
            
            if obj.dlDeltaRatio > 1
                out = obj.dlLastLambda * 1.1009;
            elseif obj.dlDeltaRatio < 1
                out = obj.dlLastLambda * 0.9084;
            else
                out = obj.dlLastLambda;
            end
            
%             out = obj.dlLastLambda * obj.dlDeltaRatio;
            
            if ~ (out < obj.dlLambdaCap && out > 0)
                
                obj.dlDeltaRatio = 1;
                obj.dlLastDelta = -1;
                
                fprintf("--->Warning! Lambda (%f) exceeded the limit (%f). Consider divergence problems.\n", out, obj.dlLambdaCap);
                fprintf("--->To avoid problems, previous lambda will be used for the next update.\n");
                fprintf("--->Consider that this problem usually happens when your model is not suitable \n");
                fprintf("--->for this optimization and you have to modify some changes. Check the hints that may help your model to perform better.");

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
           
            fprintf("\n\t->Single trial running: \n");

            try

                obj.dlUpdateParams(dlVaryList);

            catch

                fprintf("->Vary list for dsSimulate ignored due to being invalid.\n");

            end

            obj.dlSimulate();
            obj.dlCalculateOutputs(dlOutputParameters);
            fprintf("\n\t-->Simulation outputs: ");
            disp(obj.dlLastOutputs);
            
        end
        
        function dlTrain(obj, dlInputParameters, dlOutputParameters, dlTargetParameters, dlTrainOptions) 
            
            fprintf("-> Training started:\n");
            
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
                
                dlEpochs = 20;
                fprintf("-->Number of epochs was not determined in options map. Default dlEpochs = 20\n");
                
            end
            
            try
               
                dlBatchs = dlTrainOptions('dlBatchs');

            catch
                
                dlBatchs = size(dlInputParameters, 2);
                fprintf("-->Batchs was not determined in options map, default dlBatchs = %d\n", size(dlInputParameters, 2));
                
            end
            
            try

                fprintf("-->Included variable(s) on fitting: \n");
                fprintf("---> %s\n", dlTrainOptions('dlTrainIncludeList'));

            catch

                dlTrainOptions('dlTrainIncludeList') = ["_netcon", "_gAMPA", "_tau", "_gleak", "_gGABA", "_gNa", "_gK", "_gCOM", "_Eleak", "noise", "Iapp"];
                fprintf("-->Warning! No variables were included in the fitting include list.\n--->Some dynasim variables are by default assigned in:\n");
                fprintf("---> %s\n", dlTrainOptions('dlTrainIncludeList'));
            end

            try

                fprintf("-->Excluded variable(s) on fitting: \n");
                fprintf("---> %s\n", string(dlTrainOptions('dlTrainExcludeList')));

            catch

                dlTrainOptions('dlTrainExcludeList') = [];
                fprintf("-->No variables were manually excluded in the fitting exclude list.\n--->No dynasim variables are excluded by default.\n");
                
            end
            
            if ~isempty(dlTrainOptions('dlTrainIncludeList')) && ~isempty(dlTrainOptions('dlTrainExcludeList'))

                [valIX,posIX] = intersect(dlTrainOptions('dlTrainIncludeList'), dlTrainOptions('dlTrainExcludeList'));
                
                if ~isempty(posIX)

                    fprintf("-->Warning! There are %d conflicts between included and excluded lists.\n--->Following variables will be excluded despite being included: \n", length(posIX));
                    fprintf("---> %s\n",valIX);

                end
                
            end

            try
               
                dlLambda = dlTrainOptions('dlLambda');
                obj.dlLastLambda = dlLambda;        
                
            catch
                
                dlLambda = 1e-4;
                obj.dlLastLambda = dlLambda;
                fprintf("-->Lambda was not determined in options map, default dlLambda = 1e-4\n");
                
            end
            
            try
               
                dlLearningRule = dlTrainOptions('dlLearningRule');
                
            catch
                
                dlLearningRule = 'GeneralEnhancedDeltaRule';
                fprintf("-->Learning rule was not determined in options map, default dlLearningRule = 'GeneralEnhancedDeltaRule'\n");
                
            end
            
            try
               
                dlUpdateMode = dlTrainOptions('dlUpdateMode');
                
            catch
                
                dlUpdateMode = 'trial';
                fprintf("-->Update mode was not determined in options map, default dlUpdateMode = 'trial'\n");
                
            end
            
            try
               
                dlOutputLogFlag = dlTrainOptions('dlOutputLogFlag');
                
            catch
                
                dlOutputLogFlag = 0;
                
            end
            
            try
               
                dlCustomLogFlag = dlTrainOptions('dlCustomLog');
                dlCustomLogArgs = dlTrainOptions('dlCustomLogArgs');
                fprintf("-->Custom log functions will be called on every trial.\n");
                
            catch
                
                dlCustomLogFlag = "n";
                fprintf("-->No custom logs specified; Ignore this if you did not have custom logs.\n")
                
            end
            
            try
               
                dlCheckpoint = dlTrainOptions('dlCheckpoint');
                
            catch
                
                dlCheckpoint = 'true';
                fprintf("-->Checkpoint flag was not determined in options map, default dlCheckpoint = 'true'\n");
                
            end
            
            try
               
                obj.dlExcludeDiverge = dlTrainOptions('dlExcludeDiverge');
                
            catch
                
                obj.dlExcludeDiverge = 0;
                fprintf("-->ExcludeDiverge was not determined in options map, default dlExcludeDiverge = 1\n");
                
            end

            try
               
                dlUpdateVerbose = dlTrainOptions('dlUpdateVerbose');

                if dlUpdateVerbose

                    fprintf("-->dlUpdateVerbose is %d (On)\n", dlTrainOptions('dlUpdateVerbose'));

                else

                    fprintf("-->dlUpdateVerbose is %d (Off)\n", dlTrainOptions('dlUpdateVerbose'));

                end
                
            catch
                
                dlTrainOptions('dlUpdateVerbose') = 0;
                fprintf("-->dlUpdateVerbose was not determined in options map, default dlUpdateVerbose = %d (Off)\n", dlTrainOptions('dlUpdateVerbose'));
                
            end

            try
               
                dlCheckpointCoefficient = dlTrainOptions('dlCheckpointCoefficient');
                
            catch
                
                dlCheckpointCoefficient = 2.0;
                fprintf("-->Checkpoint Coefficient for optimal state saving and loading was not determined in options map, default dlCheckpointCoefficient = 2\n");
                
            end

            try
               
                dlCheckpointLengthCap = dlTrainOptions('dlCheckpointLengthCap');
                
            catch
                
                dlCheckpointLengthCap = 20;
                fprintf("-->Checkpoint LengthCap for optimal state saving and loading was not determined in options map, default value = 20\n");
                
            end
            
            if dlSimulationFlag ~= 1
               
                fprintf("-->Simulation has been manually deactivated for this run.\n");
                if dlOfflineOutputGenerator == 1
               
                    fprintf("--->Offline output generator (random) is activated. Outputs are only for debugging approaches.\n");
                
                end
                
            else
               
                if dlOfflineOutputGenerator == 1
               
                    fprintf("-->Offline output generator (random) is activated but ignored as simulation is active. \n");
                
                end
                
            end
            
            if dlAdaptiveLambda == 1
               
                fprintf("-->Adaptive lambda is active. Lambda (learning rate) will be changed based on volatility of model error.\n");
                try
                    obj.dlLambdaCap = dlTrainOptions('dlLambdaCap');
                catch
                    fprintf("--->Reminder: you have choosen adaptive lambda but forgot to determine a lambda cap. default dlLambdaCap = %f\n", obj.dlLambdaCap);
                end
            end
            
            if dlOutputLogFlag == 1
               
                fprintf("-->Outputs log will be saved.\n");
                
            end            
                
            i = 0;
            dlCurrentCheckpointLength = 0;

            while i < dlEpochs*dlBatchs
                
                fprintf("   ->Total valid trials %d (Total iterations for this model : %d)\n", max(size(obj.dlLastErrorsLog)), obj.dlTrialNumber);

                for j = 1:dlBatchs

                    dlCurrentCheckpointLength = dlCurrentCheckpointLength + 1;
                    
                    if i >= dlEpochs*dlBatchs

                        break;

                    end

                    fprintf("\t--> Trial %d of Cap. %d\n", i, dlEpochs*dlBatchs);
                    fprintf("\t--> Checkpoint. %d of Cap. %d\n", dlCurrentCheckpointLength, dlCheckpointLengthCap);
                    fprintf("\t--> Batch no. %d of %d\t\n", j, dlBatchs);
                    set(obj, 'dlTrialNumber', obj.dlTrialNumber + 1);

                    if ~isempty(dlInputParameters)

                        obj.dlUpdateParams(dlInputParameters{j});

                    end

                    if dlSimulationFlag == 1

                        obj.dlSimulate();

                        if ~isempty(dlOutputParameters)
                            obj.dlCalculateOutputs(dlOutputParameters);
                        end

                    elseif dlOfflineOutputGenerator == 1

                        obj.dlOutputGenerator();

                    end
                    
                    if ~strcmpi(dlCustomLogFlag, "n")

                        if isempty(obj.dlCustomLog)

                            obj.dlCustomLog = cell(max(size(dlCustomLogFlag)), 1);
                            obj.dlCustomLogLabel = cell(max(size(dlCustomLogFlag)), 1);

                        end

                        k = size(obj.dlCustomLog, 2);

                        for customLogCount = 1:max(size(dlCustomLogFlag))

                            logfuncname = dlCustomLogFlag(customLogCount);           
                            dlLogFuncBridge(logfuncname);
                            [cLog, cLabel] = dlLogTempFunc(obj, dlCustomLogArgs(customLogCount));

                            obj.dlCustomLog{customLogCount, k+1} = cLog;
                            obj.dlCustomLogLabel{customLogCount} = cLabel;

                        end

                    end

                    if ~isempty(dlTargetParameters)
                        obj.dlCalculateError(dlTargetParameters{j});
                        fprintf("\n\t--->Error of this trial = %f\n", obj.dlLastError);
                    end
                    
                    if dlOutputLogFlag
                        k = max(size(obj.dlOutputLog));
                        obj.dlOutputLog{k+1} = imresize(cell2mat(obj.dlOutputs), [900, 567]);
                    end
                    
                    if strcmpi(dlUpdateMode, 'trial')
                        
                        if strcmpi(dlCheckpoint, 'true')
                        
                            if obj.dlLastError < obj.dlOptimalError

                                obj.dlOptimalError = obj.dlLastError;
                                obj.dlOptimalWeightChanges = obj.dlLastWeightChanges;
                                obj.dlOptimalWeightChangesFlag = 1;
                                obj.dlSaveOptimal();

                                obj.dlUpdateError = obj.dlLastError;
                                if dlAdaptiveLambda == 1
                                    dlLambda = obj.dlAdaptiveLambda();
                                end

                                dlCurrentCheckpointLength = 0;
                                obj.dlTrainStep(dlLearningRule, dlLambda, dlTrainOptions);

                            elseif obj.dlLastError > dlCheckpointCoefficient*obj.dlOptimalError

                                disp("CoeffcCap")
                                obj.dlLoadOptimal();
                                dlCurrentCheckpointLength = 0;

                            elseif isnan(obj.dlLastError)

                                disp("NaNCap")
                                obj.dlLoadOptimal();
                                dlCurrentCheckpointLength = 0;

                            elseif dlCheckpointLengthCap < dlCurrentCheckpointLength

                                disp("LengthCap*");
                                obj.dlLoadOptimalLX();
                                dlCurrentCheckpointLength = 0;

                            else
                                
                                obj.dlUpdateError = obj.dlLastError;
                                if dlAdaptiveLambda == 1
                                    dlLambda = obj.dlAdaptiveLambda();
                                end
                                obj.dlTrainStep(dlLearningRule, dlLambda, dlTrainOptions);

                            end
                            
                        else
                            
                            obj.dlUpdateError = obj.dlLastError;
                            if dlAdaptiveLambda == 1
                                dlLambda = obj.dlAdaptiveLambda();
                            end
                            obj.dlTrainStep(dlLearningRule, dlLambda, dlTrainOptions);
                        
                        end
                        
                    end

                    i = size(obj.dlErrorsLog, 2) - obj.dlCurrentSessionValidTrials;
                    
                end
                
                if ~isempty(obj.dlErrorsLog)
                    if numel(obj.dlErrorsLog) > 2
                        dlAvgError = mean(obj.dlErrorsLog(end-2:end));
                    else
                        dlAvgError = mean(obj.dlErrorsLog);
                    end
                    fprintf("\t-->Epoch's Average Error = %f, Last lambda = %.14f\n", dlAvgError, dlLambda);
                end
                
                if strcmpi(dlCheckpoint, 'true')
                    
                    if dlAvgError < obj.dlOptimalError

                        obj.dlOptimalError = dlAvgError;
                        obj.dlOptimalWeightChanges = obj.dlLastWeightChanges;
                        obj.dlOptimalWeightChangesFlag = 1;
                        obj.dlSaveOptimal();

                        dlCurrentCheckpointLength = 0;

                        if strcmpi(dlUpdateMode, 'batch')

                            obj.dlUpdateError = dlAvgError;

                            if dlAdaptiveLambda == 1
                                dlLambda = obj.dlAdaptiveLambda();
                            end

                            obj.dlTrainStep(dlLearningRule, dlLambda, dlTrainOptions);

                        end

                    elseif dlAvgError > dlCheckpointCoefficient*obj.dlOptimalError

                        obj.dlLoadOptimal();
                        dlCurrentCheckpointLength = 0;

                    elseif dlCheckpointLengthCap < dlCurrentCheckpointLength

                        disp("LengthCap*");
                        obj.dlLoadOptimalLX();
                        dlCurrentCheckpointLength = 0;

                    elseif isnan(obj.dlLastError)

                            disp("NaNCap")
                            obj.dlLoadOptimal();
                            dlCurrentCheckpointLength = 0;

                    else

                        if strcmpi(dlUpdateMode, 'batch')

                            obj.dlUpdateError = dlAvgError;
                            if dlAdaptiveLambda == 1
                                dlLambda = obj.dlAdaptiveLambda();
                            end
                            obj.dlTrainStep(dlLearningRule, dlLambda, dlTrainOptions);

                        end

                    end
                    
                else
                   
                    if strcmpi(dlUpdateMode, 'batch')
                        if exist('dlAvgError', 'var')
                            obj.dlUpdateError = dlAvgError;
                            if dlAdaptiveLambda == 1
                                dlLambda = obj.dlAdaptiveLambda();
                            end
                        end
                        obj.dlTrainStep(dlLearningRule, dlLambda, dlTrainOptions);
                    end
                    
                end
                
                i = size(obj.dlErrorsLog, 2) - obj.dlCurrentSessionValidTrials;

            end

            obj.dlCurrentSessionValidTrials = size(obj.dlErrorsLog, 2);
            
        end
        
        function dlTrainStep(obj, dlLearningRule, dlLambda, dlTrainOptions)

            p = load([obj.dlPath, '/params.mat']);
            lab = fieldnames(p.p);
            val = struct2cell(p.p);

            %%% Local plasticity rules

            try

                Local_LR_Params = dlTrainOptions('Local_LR_Params');
                localFlag = true;

            catch

                localFlag = false;

            end

            try

                dlEnhancedMomentum = dlTrainOptions('dlEnhancedMomentum');

            catch

                dlEnhancedMomentum = 0.1;
                disp("----->Deflt. used Momntm. = 0.1");

            end

            if localFlag

                kernel = 'L'; % 'E'; % 'G'; % Laplacian, Epanechnikov or Gaussian
                uscaling = 1e3; % unit scaling from kernel regression from kHz to Hz
                kwidth = 500; % width of kernel regression in ms
                Ts = 1;  % subsampling period in ms
                flag_interp = 1;

                for iLocalParams = 1:numel(Local_LR_Params) % e.g., {E->E, I->E}
                    LocalParams = Local_LR_Params{iLocalParams};

                    learningRules = LocalParams.learningRules;
                    connection_type = LocalParams.connection_type;
                    source = LocalParams.source;
                    target = LocalParams.target;
                    fr_norm = LocalParams.fr_norm;

                    w_min = LocalParams.w_min;
                    w_max = LocalParams.w_max;
                    voltage = LocalParams.voltage;

                    % Connection index
                    i_conn = find(ismember([obj.dlGraph.edges.target], target) & ismember([obj.dlGraph.edges.source], source) & ismember([obj.dlGraph.edges.type], connection_type));

                    % Indices into dlGraph for presynaptic and postsynaptic populations
                    i_pre = obj.dlGraph.IndexMap(obj.dlGraph.edges(i_conn).source);
                    i_post = obj.dlGraph.IndexMap(obj.dlGraph.edges(i_conn).target);

                    % Indices to V
                    i_Vpre = find(strcmpi(string(obj.dlVariables), obj.dlGraph.vertices(i_pre).name + "_" + voltage));
                    i_Vpost = find(strcmpi(string(obj.dlVariables), obj.dlGraph.vertices(i_post).name + "_" + voltage));

                    % Voltages
                    Vpre = obj.dlOutputs{i_Vpre};
                    Vpost = obj.dlOutputs{i_Vpost};

                    % retrieve iFRs:
                    t = obj.dlOutputs{1};
                    raster = computeRaster(t, Vpre);
                    x = zeros(size(Vpre));
                    if size(raster, 1) > 0
                        % pool = 1:size(Vpre, 2);
                        % x = uscaling*NWKraster(t, raster, pool, kwidth, Ts, flag_interp, kernel);
                        for ipool = 1:size(Vpre, 2)
                            x(:,ipool) = uscaling*NWKraster(t, raster, ipool, kwidth, Ts, flag_interp, kernel);
                        end
                    end
                    x = x/fr_norm;
                    X = mean(x,1);

                    raster = computeRaster(t, Vpost);
                    y = zeros(size(Vpost));
                    if size(raster, 1) > 0
                        % pool = 1:size(Vpost, 2);
                        % y = uscaling*NWKraster(t, raster, pool, kwidth, Ts, flag_interp, kernel);
                        for ipool = 1:size(Vpost, 2)
                            y(:,ipool) = uscaling*NWKraster(t, raster, ipool, kwidth, Ts, flag_interp, kernel);
                        end
                    end
                    y = y/fr_norm;
                    Y = mean(y,1);

                    % Label of netcon in params.mat
                    lkey = obj.dlGraph.vertices(i_post).name + "_" + obj.dlGraph.vertices(i_pre).name + "_" + obj.dlGraph.edges(i_conn).type + "_netcon";

                    % Connectivity matrix
                    ind_w = find(strcmpi(string(lab), lkey)); % Find index of W in params.mat labels
                    w = val{ind_w}; % previous connectivity matrix

                    if isempty(find(strcmpi(string(obj.dlWeightsVariables), lkey)))
                        obj.dlWeightsVariables{numel(obj.dlWeightsVariables)+1} = lkey;
                        obj.dlWeightsValues{numel(obj.dlWeightsValues)+1}(:,:,1) = w;
                    end

                    local_delta = zeros(size(w));

                    for iLearningRule = 1:numel(learningRules) % e.g., {ACh, NE, Hebbian} or {Inhibitory}
                        learningRule = learningRules{iLearningRule};

                        %%% What follows are 4 Local Learning Rules (3 for E-cells {Ach, NE, and Hebbian} and 1 for I-cells {Inhibitory}) from Aljadeff et al., arXiv 2019:
                        % Aljadeff et al. Cortical credit assignment by Hebbian, neuromodulatory and inhibitory plasticity. arXiv:1911.00307, 2019

                        if strcmpi(learningRule, 'ACh')
                            % ACh learning rule based on Aljadeff et al., arXiv 2019 (LTP and LTD)
                            % x: presynaptic activity
                            % y: postsynaptic activity
                            % y_ref: reference firing rate for neuromodulation plasticity (with y_ref~0, no LTD even if y = 0)
                            % rho: input-output pairing probability, e.g. based on selectivity similarity
                            % Y: binary signal that enables ACh neuromodulation as a whole
                            % eta: if pairing is enabled, eta activates ACh plasticity based on rho
                            % alpha: learning rate for ACh plasticity
                            % beta: LTD/LTP scaling factor
                            % single expression for balanced LTD/LTP ('*' represents matrix multiplication):
                            % delta = eta·alpha·(y-y_ref)*(x - beta·f/(1-f)·(1-x))

                            if ~LocalParams.Y_ACh
                                continue
                            end

                            rho_ACh = LocalParams.rho_ACh;
                            alpha_ACh = LocalParams.alpha_ACh;
                            y_ref_ACh = LocalParams.y_ref_ACh;
                            beta_ACh = LocalParams.beta_ACh;
                            f_ref = LocalParams.f_ref;

                            % eta_ACh = rand(1, size(w,2)) <= rho_ACh;
                            eta_ACh_name = LocalParams.eta_ACh;
                            i_eta_ACh = find(strcmpi(string(obj.dlVariables), obj.dlGraph.vertices(i_post).name + "_" + eta_ACh_name));
                            eta_ACh = obj.dlOutputs{i_eta_ACh};

                            local_delta = local_delta + alpha_ACh*(X' - beta_ACh*f_ref/(1-f_ref).*(1-X'))*(eta_ACh.*(Y-y_ref_ACh));

                        elseif strcmpi(learningRule, 'NE')
                            % NE learning rule based on Aljadeff et al., arXiv 2019 (only LTP and not stimulus specific)
                            % x: presynaptic activity
                            % y: postsynaptic activity
                            % y_ref: reference firing rate for neuromodulation plasticity
                            % rho: input-output pairing probability, e.g. based on selectivity similarity
                            % eta: enables NE plasticity based on rho
                            % alpha: learning rate for NE plasticity
                            % expression:
                            % delta = eta·alpha·(y-y_ref)*x

                            rho_NE = LocalParams.rho_NE;
                            alpha_NE = LocalParams.alpha_NE;
                            y_ref_NE = LocalParams.y_ref_NE;

                            % eta_NE = rand(1, size(w,2)) <= rho_NE;
                            eta_NE_name = LocalParams.eta_NE;
                            i_eta_NE = find(strcmpi(string(obj.dlVariables), obj.dlGraph.vertices(i_post).name + "_" + eta_NE_name));
                            eta_NE = obj.dlOutputs{i_eta_NE};

                            local_delta = local_delta + alpha_NE*X'*(eta_NE.*(Y-y_ref_NE));

                        elseif strcmpi(learningRule, 'Hebbian')
                            % Hebbian learning rule based on Aljadeff et al., arXiv 2019 (not gated)
                            % x: presynaptic activity
                            % y: postsynaptic activity
                            % f_ref: reference firing rate for neuromodulation plasticity
                            % alpha: learning rate for Hebbian plasticity
                            % expression:
                            % delta = alpha·(y-f_ref)*x

                            alpha_Hebbian = LocalParams.alpha_Hebbian;
                            y_ref_Hebbian = LocalParams.y_ref_Hebbian;

                            local_delta = local_delta + alpha_Hebbian*X'*(Y-y_ref_Hebbian);

                        elseif strcmpi(learningRule, 'Inhibitory')
                            % learning rule for I-cells based on Aljadeff et al., arXiv 2019 (based on detailed E-I balance)
                            % x: presynaptic activity
                            % excitatory current iE: all excitatory currents to E-cells
                            % inhibitory current iI: all inhibitory currents to E-cells
                            % E/I balance line (reference): iI = a·iE + b (with a <= 1)
                            % alpha: learning rate for Inhibitory plasticity
                            % expression:
                            % delta = alpha·((a·iE + b)-iI)*x

                            a_Inhibitory = LocalParams.a_Inhibitory;
                            b_Inhibitory = LocalParams.b_Inhibitory;
                            alpha_Inhibitory = LocalParams.alpha_Inhibitory;

                            % Variables
                            iE_name = LocalParams.iE;
                            iI_name = LocalParams.iI;

                            % Indices to postsynaptic iE and iI
                            i_iE = find(strcmpi(string(obj.dlVariables), [obj.dlGraph.vertices(i_post).name, '_', iE_name]));
                            i_iI = find(strcmpi(string(obj.dlVariables), [obj.dlGraph.vertices(i_post).name, '_', iI_name]));

                            % Currents
                            iE = mean(obj.dlOutputs{i_iE},1);
                            iI = mean(obj.dlOutputs{i_iI},1);

                            local_delta = local_delta + alpha_Inhibitory*X'*((a_Inhibitory*iE + b_Inhibitory) - iI);
                        else
                            error('unknown local plasticity rule');
                        end
                    end
                    % Update weights
                    w = w + local_delta;
                    if any(w < w_min)
                        warning('weights below lower bound, rectifying')
                        w(w < w_min) = w_min;
                    end
                    if any(w > w_max)
                        warning('weights above upper bound, rectifying')
                        w(w > w_max) = w_max;
                    end
                    val{ind_w} = w;
                    i_w = find(strcmpi(string(obj.dlWeightsVariables), lkey));
                    obj.dlWeightsValues{i_w}(:,:,end+1) = w;
               end
           end

           %%% Non-local plasticity rules

           error = obj.dlUpdateError;
           obj.dlLastLambda = dlLambda;
           l = find(contains(lab, dlTrainOptions('dlTrainIncludeList'))); 
           lg = find(contains(lab, dlTrainOptions('dlTrainIncludeList')));

           deltaL = 0;

           if strcmpi(dlLearningRule, 'DeltaRule')

               for i = l'

                   rng('shuffle');
                   w = val{i, 1};
                   delta = (randn(size(w)))*error*dlLambda;

                   if ~contains(lab{i, 1}, 'IO')

                       wn = w + delta;

                       wn(wn < .1) = .1;
                       wn(wn > .9) = .9;
                       val{i, 1} = wn;
                       deltaL = deltaL + sum(sum(abs(delta)));

                   end

               end
                
           elseif strcmpi(dlLearningRule, 'BioDeltaRule')
            
               restrictedCoefCnt = 0;

               for i = l'

                   rng('shuffle');
                   w = val{i, 1};
                   delta = (1-w*0.9).*(rand(size(w))-0.5)*error*dlLambda;
                   obj.dlLastWeightChanges{i} = delta;

                   try

                        excludeList = contains(lab, dlTrainOptions('dlTrainExcludeList'));

                   catch

                        excludeList = contains(lab, '');
                        excludeList = excludeList * 0;

                   end

                   try

                        restrictedList = contains(lab, dlTrainOptions('dlTrainRestrictList'));
                        restrictedCoef = dlTrainOptions('dlTrainRestrictCoef');

                   catch

                        restrictedList = contains(lab, '');
                        restrictedList = restrictedList * 0;

                   end

                   if ~excludeList(i)

                    if restrictedList(i)

                        restrictedCoefCnt = restrictedCoefCnt + 1;
                        wn = w - delta*restrictedCoef{restrictedCoefCnt};
                    
                    else

                        wn = w - delta;

                    end
                    
                    wn(wn < 0.1) = 0.1;
                    wn(wn > .94) = .94;
                    val{i, 1} = wn;
                    deltaL = deltaL + sum(sum(abs(delta)));

                   end

               end
           
           elseif strcmpi(dlLearningRule, 'EnhancedDeltaRule')

            restrictedCoefCnt = 0;
            for i = l'

                rng('shuffle');
                w = val{i, 1};
                delta = (1-w*0.9).*(rand(size(w))-0.5)*error*dlLambda;
                obj.dlLastWeightChanges{i} = delta;
                
                if obj.dlOptimalWeightChangesFlag
                
                   try
                
                        delta = delta*(1-dlEnhancedMomentum) + obj.dlOptimalWeightChanges{i}*dlEnhancedMomentum;
                 
                   catch
                
                       if i == l(1)
                           fprintf("<q107>");
                       end

                   end
                
                end
                
                try
                
                    excludeList = contains(lab, dlTrainOptions('dlTrainExcludeList'));
                
                catch
                
                    excludeList = contains(lab, '');
                    excludeList = excludeList * 0;
                
                end
                
                try
                
                    restrictedList = contains(lab, dlTrainOptions('dlTrainRestrictList'));
                    restrictedCoef = dlTrainOptions('dlTrainRestrictCoef');
                
                catch
                
                    restrictedList = contains(lab, '');
                    restrictedList = restrictedList * 0;
                
                end

                    if ~excludeList(i)
                    
                        if restrictedList(i)
                        
                            restrictedCoefCnt = restrictedCoefCnt + 1;
                            wn = w - delta*restrictedCoef{restrictedCoefCnt};
                        
                        else
                        
                            wn = w - delta;
                        
                        end
                        
                           wn(wn < 0.1) = 0.1;
                           wn(wn > .94) = .94;
                           val{i, 1} = wn;
                           deltaL = deltaL + sum(sum(abs(delta)));
                        
                    end

              end

            elseif strcmpi(dlLearningRule, 'GeneralEnhancedDeltaRule')

                restrictedCoefCnt = 0;
    
                for i = lg'
    
                    rng('shuffle');
                    w = val{i, 1};

                    if contains(lab{i, 1}, '_netcon')

                        delta = (1-w*0.9).*(rand(size(w))-0.5)*error*dlLambda;

                    else

                        delta = w*(rand(size(w))-0.5)*error*dlLambda;

                    end

                    obj.dlLastWeightChanges{i} = delta;
                    
                    if obj.dlOptimalWeightChangesFlag
                    
                       try
                    
                            delta = delta*(1-dlEnhancedMomentum) + obj.dlOptimalWeightChanges{i}*dlEnhancedMomentum;
                     
                       catch

                           if i == lg(1)
                               fprintf("\n<q107:Not enough differential initial values for momentum. \n-> ..." + ...
                               "Establishing initial values takes several trials; This is not an error.>\n");
                           end

                       end
                    
                    end
                    
                    try
                    
                        excludeList = contains(lab, dlTrainOptions('dlTrainExcludeList'));
                    
                    catch
                    
                        excludeList = contains(lab, '');
                        excludeList = excludeList * 0;
                    
                    end
                    
                    try
                    
                        restrictedList = contains(lab, dlTrainOptions('dlTrainRestrictList'));
                        restrictedCoef = dlTrainOptions('dlTrainRestrictCoef');
                    
                    catch
                    
                        restrictedList = contains(lab, '');
                        restrictedList = restrictedList * 0;
                    
                    end
    
                    if ~excludeList(i)
                    
                        if restrictedList(i)
                        
                            restrictedCoefCnt = restrictedCoefCnt + 1;
                            wn = w - delta*restrictedCoef{restrictedCoefCnt};
                        
                        else
                        
                            wn = w - delta;
                        
                        end
            
                        if contains(lab{i, 1}, "_netcon")

                            wn(wn < 0.1) = 0.1;
                            wn(wn > .94) = .94;

                        elseif contains(lab{i, 1}, "_E")

                            wn(abs(wn) > 1e3) = sign(wn)*9e2; % Too large elimination for stability
                        
                        else

                            wn(abs(wn) < 1e-3) = 9e-2; % Too small or negative elimination for stability
                            wn(abs(wn) > 1e+4) = 9e+3*sign(wn); % Too high values elimination for stability

                        end

                        val{i, 1} = wn;

                        if dlTrainOptions('dlUpdateVerbose')

                            fprintf("\n -----> Updated: %s : %f \n", lab{i, 1}, val{i, 1});

                        end

                        deltaL = deltaL + sum(sum(abs(delta)));
                        
                    end
    
                end

           elseif strcmpi(dlLearningRule, 'UncertaintyReduction')

               error = obj.dlLastError; %TODO ask Hamed about the difference wrt obj.dlUpdateError

               metaLR0 = obj.dlMetaLR;
               if isscalar(metaLR0)
                   obj.error0 = error;
               end
               uncBaseline = 0.5; % reference value for uncertainty reduction (point of maximum uncertainty)
               scalingFactor = max([uncBaseline, 1-uncBaseline]); % used to keep uncReduct in [0,1]
               stochasticFactor = 0.2; % stochastic modulation

               for i = l'
                   rng('shuffle'); %% TODO we shouldn't shuffle all the time

                   w = val{i};

                   % lambda update from previous w
                   if isscalar(metaLR0) % first time
                       alpha = metaLR0*ones(size(w));
                   else % subsequent times
                       alpha = metaLR0{i};
                   end
                   uncReduct = abs(w-uncBaseline)/scalingFactor; % Uncertainty reduction in [0,1]
                   % adapting lambda based on uncertainty reduction (the lower the uncertainty, the faster it adapts)
                   alpha = alpha + obj.dlMetaMu*(uncReduct-alpha); % adapting lambda based on uncertainty reduction

                   % w update based on the new lambda
                   delta = (1 + stochasticFactor*randn(size(w))).*alpha.*((1-exp(-abs(error/obj.error0))) - w); % stochastic delta
                   w = w + delta;
                   % rectifying values that are out of the [0,1] bounds
                   if any(w < 0)
                       warning('weights below lower bound, rectifying')
                       w(w < 0) = 0;
                   end
                   if any(w > 1)
                       warning('weights above upper bound, rectifying')
                       w(w > 1) = 1;
                   end

                   % saving
                   val{i} = w;
                   metaLR{i} = alpha;

                   deltaL = deltaL + sum(abs(delta(:)));
               end
               obj.dlMetaLR = metaLR;

           elseif strcmpi(dlLearningRule, 'RWDeltaRule')
            
               disp("TODO Rascorla-Wagner delta rule");
                
           elseif strcmpi(dlLearningRule, 'NewRule')
            
                % This section is only implemented to show how to use the
                % obj.dlGraph for local learning purposes. >>>>>>>>>>>>>>
                connectionCount = size(obj.dlGraph.edges, 2); % Number of edges or connections in model
                populationCount = size(obj.dlGraph.vertices, 2); % Number of vertices or populations in model
                
                fprintf("\n->NewRule for learning method is only an example of dlGraph and does nothing meaningfull\n");

                for i = 1:connectionCount % Traverse all connections (edges)
                    
                    if contains(obj.dlGraph.edges(i).type, 'GABA')
                        fprintf("\n %s to %s is Inhibitory Connection!", obj.dlGraph.edges(i).source, obj.dlGraph.edges(i).target); % Displays if the connection type is GABA/inhibitory
                    elseif contains(obj.dlGraph.edges(i).type, 'AMPA')
                        fprintf("\n %s to %s is Excitatory Connection!", obj.dlGraph.edges(i).source, obj.dlGraph.edges(i).target); % Displays if the connection type is AMPA/excitatory
                    end
                    
                    v = obj.dlGraph.IndexMap(obj.dlGraph.edges(i).source); % find source ...
                    u = obj.dlGraph.IndexMap(obj.dlGraph.edges(i).target); % find connection(i) target's index in vertices
                    
                    IndexInOutputs = find(strcmpi(string(obj.dlVariables), obj.dlGraph.vertices(u).name + "_V")); % find its voltage in the current outputs
                    lkey = obj.dlGraph.vertices(u).name + "_" + obj.dlGraph.vertices(v).name + "_" + obj.dlGraph.edges(i).type + "_netcon"; % netcon's label in params.mat
                    WIndexInParams = find(strcmpi(string(lab), lkey)); % Find index of W in params.mat labels
                    
                    V = obj.dlOutputs{1, IndexInOutputs}; % extract its V
                    v = mean(mean(V));
                    delta = (1-v).*(randn(size(v)))*error*dlLambda; % Do some calculation
                    
                    try
                        w = val{WIndexInParams, 1};
                        wn = w + delta;
                        wn(wn < 0) = 0; % Apply some treshold (just for example)
                        wn(wn > 1) = 1;
                        val{WIndexInParams, 1} = wn; % Put new W in the old W
                    catch
%                         fprintf("\nNo specific type for connection, probably input or output"); % Commented for shorter log
                    end
                        
                end
                
                for i = 1:populationCount % Traverse all populations (vertices)
                    
                    rng('shuffle');
                    
                    preSynaptic = obj.dlGraph.vertices(i).PreIndices; % find vertex(i) presynaptic indices in vertices
                    postSynaptic = obj.dlGraph.vertices(i).PostIndices; % find vertex(i) postsynaptic indices in vertices
                    preSynapticConnections = obj.dlGraph.vertices(i).PreEdgeIndices; % find vertex(i) presynaptic edge (connection) indices in obj.dlGraph.edges
                    postSynapticConnections = obj.dlGraph.vertices(i).PostEdgeIndices; % find vertex(i) postsynaptic edge (connection) indices in obj.dlGraph.edges
                    
                    PreIndexInOutputs = [];
                    PostIndexInOutputs = [];
                    
                    for j = preSynaptic
                        PreIndexInOutputs = [PreIndexInOutputs, find(strcmpi(string(obj.dlVariables), obj.dlGraph.vertices(j).name + "_V"))]; % find presynaptic voltages' indices in the current outputs
                    end
                    
                    for j = postSynaptic
                        PostIndexInOutputs = [PostIndexInOutputs, find(strcmpi(string(obj.dlVariables), obj.dlGraph.vertices(j).name + "_V"))]; % find postsynaptic voltages' indices in the current outputs
                    end
                    
                    w = 0;
                    for j = PreIndexInOutputs
                        w = w + mean(mean(val{j, 1})); % use and do some operation on Presynaptic Ws
                    end
                    for j = PostIndexInOutputs
                        w = w + mean(mean(val{j, 1})); % use and do some operation on Presynaptic Ws
                    end
                    
                    delta = (1-w).*(randn(size(w)))*error*dlLambda; % Do some calculation
                    
                    if ~contains(lab{IndexInOutputs, 1}, 'IO') % Check if it has/has'nt some keyword
                        wn = w + delta;
                    end
                    
                    wn(wn < 0) = 0; % Apply some treshold (just for example)
                    wn(wn > 1) = 1;
                    val{IndexInOutputs, 1} = wn; % Put new W in the old W
                    
                end

            elseif strcmpi(dlLearningRule, 'none')

                % disp('not using a dlLearningRule, most likely using a local plasticity rule instead')

            else
                
                disp("TODO train step and learning 'else' part");
                
            end
            
            if obj.dlLastDelta < 0
            
                obj.dlLastDelta = deltaL;
                    
            else
                
                dlErrorChangesPenalty = ((obj.dlLastError - obj.dlOptimalError) / obj.dlOptimalError)^0.5;
                obj.dlDeltaRatio = max(min((obj.dlLastDelta / deltaL)^0.5, 2), dlErrorChangesPenalty);
                obj.dlLastDelta = deltaL;
                
            end
            
            q = cell2struct(val, lab);
            % disp(val);
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
            x = zeros(1, n);
            xm = zeros(1, n);
            
            for i = 0:n-dlBatchs
                x(i+1) = mean(obj.dlErrorsLog(i+1:i+dlBatchs));
                xm(i+1) = min(obj.dlErrorsLog(i+1:i+dlBatchs));
            end
            
            % plot(x, "DisplayName", "Mean error", "LineStyle", "--");
            % grid("on");hold("on");
            % plot(xm, "DisplayName", "Min error");
            % xlabel("Trials");
            % ylabel("Log10(Loss)");
            % legend;
            % title("Errors in batchs");

            plot(log10(x), "DisplayName", "Mean error", "LineStyle", "--");
            grid("on");hold("on");
            plot(log10(xm), "DisplayName", "Min error");
            xlabel("Trials");
            ylabel("Log10(Loss)");
            legend;
            title("Errors in batchs");
            
        end
        
        function dlUpdateParams(obj, map) 
            
            fprintf("Updating parameters of %s > \n", obj.dlPath);
            dsParamsModifier('dlTempFuncParamsChanger.m', map);
            dlTempFuncParamsChanger(obj.dlPath);
            
        end
        
        function dlResetTraining(obj)
           
            obj.dlTrialNumber = 0;
            obj.dlOptimalError = 1e9;
            obj.dlOutputLog = [];
            obj.dlErrorsLog = [];
            
            obj.dlLastErrorsLog = [];
            obj.dlCustomLog = [];
            obj.dlCustomLogLabel = [];
            obj.dlCurrentSessionValidTrials = 0;

            obj.dlLastOutputLog = {};
            
        end
        
        function dlSaveOptimal(obj)
           
            obj.dlLastErrorsLog = obj.dlErrorsLog;
            obj.dlLastCustomLog = obj.dlCustomLog;
            obj.dlLastOutputLog = obj.dlOutputLog;
            obj.dlOutputs = [];

            obj.dlSaveCheckPoint('/Optimal');
            
        end
        
        function dlLoadOptimal(obj)
            
            obj.dlLoadCheckPoint('/Optimal');
 
        end

        function dlLoadOptimalLX(obj)

    
            obj.dlLoadCheckPointLX('/Optimal');
        
        end

        function dlLoadOptimalLC(obj)

            try

                obj.dlLoadCheckPoint('/Optimal');

            catch

                fprintf("--->No oprimal file exists. first run a training session with an active checkpoint flag to save an optimal checkpoint.\n");
          
            end

        end
        
        function dlGraphConstructor(obj) 
            
            fprintf("--->Constructing graph ...");
            e = size(obj.dlModel.connections, 2);
            v = size(obj.dlModel.populations, 2);
            
            obj.dlGraph.edges = {};
            obj.dlGraph.vertices = {};
            obj.dlGraph.IndexMap = containers.Map();
            
            for i = 1:v
               
                obj.dlGraph.IndexMap(obj.dlModel.populations(i).name) = i;
                obj.dlGraph.vertices(i).name = obj.dlModel.populations(i).name;
%               IndexInOutputs = find(strcmpi(string(obj.dlVariables), obj.dlGraph.vertices(i).name + "_V"));
                obj.dlGraph.vertices(i).PreIndices = 0;
                obj.dlGraph.vertices(i).PostIndices = 0;
                
            end
            
            for i = 1:e
                
                try
                    
                    s = string(split(obj.dlModel.connections(i).direction, "->"));
                    obj.dlGraph.edges(i).source = s(1);
                    obj.dlGraph.edges(i).target = s(2);
                
                catch
                    
                    obj.dlGraph.edges(i).source = obj.dlModel.connections(i).source;
                    obj.dlGraph.edges(i).target = obj.dlModel.connections(i).target;
                
                end
                
                if ~obj.dlGraph.vertices(obj.dlGraph.IndexMap(obj.dlGraph.edges(i).source)).PostIndices
                    obj.dlGraph.vertices(obj.dlGraph.IndexMap(obj.dlGraph.edges(i).source)).PostIndices = obj.dlGraph.IndexMap(obj.dlGraph.edges(i).target);
                    obj.dlGraph.vertices(obj.dlGraph.IndexMap(obj.dlGraph.edges(i).source)).PostEdgeIndices = i;
                else
                    obj.dlGraph.vertices(obj.dlGraph.IndexMap(obj.dlGraph.edges(i).source)).PostIndices = cat(2, obj.dlGraph.vertices(obj.dlGraph.IndexMap(obj.dlGraph.edges(i).source)).PostIndices, obj.dlGraph.IndexMap(obj.dlGraph.edges(i).target));
                    obj.dlGraph.vertices(obj.dlGraph.IndexMap(obj.dlGraph.edges(i).source)).PostEdgeIndices = cat(2, obj.dlGraph.vertices(obj.dlGraph.IndexMap(obj.dlGraph.edges(i).source)).PostEdgeIndices, i);
                end

                if ~obj.dlGraph.vertices(obj.dlGraph.IndexMap(obj.dlGraph.edges(i).target)).PreIndices
                    obj.dlGraph.vertices(obj.dlGraph.IndexMap(obj.dlGraph.edges(i).target)).PreIndices = obj.dlGraph.IndexMap(obj.dlGraph.edges(i).source);
                    obj.dlGraph.vertices(obj.dlGraph.IndexMap(obj.dlGraph.edges(i).target)).PreEdgeIndices = i;
                else
                    obj.dlGraph.vertices(obj.dlGraph.IndexMap(obj.dlGraph.edges(i).target)).PreIndices = cat(2, obj.dlGraph.vertices(obj.dlGraph.IndexMap(obj.dlGraph.edges(i).target)).PreIndices, obj.dlGraph.IndexMap(obj.dlGraph.edges(i).source));
                    obj.dlGraph.vertices(obj.dlGraph.IndexMap(obj.dlGraph.edges(i).target)).PreEdgeIndices = cat(2, obj.dlGraph.vertices(obj.dlGraph.IndexMap(obj.dlGraph.edges(i).target)).PreEdgeIndices, i);
                end
                
                try
                    
                    obj.dlGraph.edges(i).type = obj.dlModel.connections(i).mechanism_list;
                
                catch
                    
                    obj.dlGraph.edges(i).type = "null";
                
                end
                
            end
            
            fprintf(" Done.\n");
            
        end
        
    end
end
