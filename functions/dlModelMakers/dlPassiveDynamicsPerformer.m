function m = dlPassiveDynamicsPerformer(RemakeFlag, ResetOptimalError, ModelName, ModelParameters, model_size_id, performance_coefficient, tune_flag, epochs)
    
    %%% Call Laminar Cortex Constructor Functions
    dsCellLaminar = dlLaminarCortexNetNL(ModelParameters, ModelName); % Laminar PFC model with specific parameters
    
    %%% Finalization
    dsModel = dsCellLaminar;
    
    % Create DynaLearn Class (remake)
    if RemakeFlag
        fprintf("\n->Remake flag is on. Generating new model (will replace the previous model if names are same)");
        m = DynaLearn(dsModel, char("dlModels/dlPredictiveCorticalCircuitModelNL" + string(model_size_id)), 'mex', ModelName); % ~10 min or less, MEXGEN or < 20 sec, RAWGEN.
        m.dlSave(); % < 1sec
    end

    % Load DynaLearn Class
    m = DynaLearn(); % ~ 1sec

    try
        m = m.dlLoad(char("dlModels/dlPredictiveCorticalCircuitModelNL" + string(model_size_id))); % ~ 10sec, New larger model; keeping track of its activity in Gamma/Beta **
    catch
        fprintf("\n->Failed to load object. If this object does not exist, try RemakeFlag=1.");
        fprintf("\n-->Otherwise, there is a problem with loading object or access to its repository.");
    end

    ModelName = char("_" + ModelName);
    [trialParams1, trialParams2, trialParams3] = dlNullPattern(ModelName);
    
    outputParams = [{['DeepE', ModelName, '_V'], 1:floor(ModelParameters.NeDeep/3), [200 400] ...
        , 'afr'}; {['DeepE', ModelName, '_V'],ceil(ModelParameters.NeDeep/3):floor(2*ModelParameters.NeDeep/3), ...
        [200 400], 'afr'}; {['DeepE', ModelName, '_V'], ceil(2*ModelParameters.NeDeep/3):ModelParameters.NeDeep, [200 400], 'afr'}; ...
        {['supE', ModelName, '_V'], 1:ModelParameters.NeSuperficial, [50 700], 'afr'}; ...
        {['MidE', ModelName, '_V'], 1:ModelParameters.NeMid, [50 700], 'afr'}; ...
        {['DeepE', ModelName, '_V'], 1:ModelParameters.NeDeep, [50 700], 'afr'}; ...
        {['supIPV', ModelName, '_V'], 1:ModelParameters.NPvSuperficial, [50 700], 'afr'}; ...
        {['MidIPV', ModelName, '_V'], 1:ModelParameters.NPvMid, [50 700], 'afr'}; ...
        {['DeepIPV', ModelName, '_V'], 1:ModelParameters.NPvDeep, [50 700], 'afr'}];
    
    perf_c = performance_coefficient;
    targetParams1 = [{'RPenalty', 4:9, 10, 1.01, 8, 32, 40, 120}; {'Compare', [1, 2], 0, perf_c, 0, 0, 0, 0}; {'Compare', [1, 3], 0, perf_c, 0, 0, 0, 0}]; % A 
    targetParams2 = [{'RPenalty', 4:9, 10, 1.01, 8, 32, 40, 120}; {'Compare', [2, 1], 0, perf_c, 0, 0, 0, 0}; {'Compare', [2, 3], 0, perf_c, 0, 0, 0, 0}]; % B
    targetParams3 = [{'RPenalty', 4:9, 10, 1.01, 8, 32, 40, 120}; {'Compare', [3, 1], 0, perf_c, 0, 0, 0, 0}; {'Compare', [3, 2], 0, perf_c, 0, 0, 0, 0}]; % C
    
    dlInputParameters = {trialParams1, trialParams2, trialParams3};
    dlTargetParameters = {targetParams1, targetParams2, targetParams3};
    dlOutputParameters = outputParams;
    
    TBdata = dlTrialBlockGenerator(dlInputParameters, dlTargetParameters, epochs, epochs);
    
    dlTrainOptions = containers.Map(); % Train options; MUST be a map data structure
    dlTrainOptions('dlEpochs') = 100; % % Number of epochs (A.K.A total iterations)
    dlTrainOptions('dlBatchs') = 3; % If a scenario requires the training to be based on a group parameter (e.g mean of errors) use a dlBatch > 1 and set update mode later to batch. 
    dlTrainOptions('dlLambda') = 5e-2; % Higher lambda means more changes based on error, lower may cause model to learn slower or nothing.
        
    dlTrainOptions('dlCheckpoint') = 'true'; % If current step's error is higher based on a threshold, reload last optimal state and continue from that point
    dlTrainOptions('dlCheckpointCoefficient') = 2.047; % A.K.A exploration rate
    dlTrainOptions('dlCheckpointLengthCap') = 7; % If more than 7 steps with no progress passed, return to last checkpoint.
    dlTrainOptions('dlUpdateMode') = 'batch'; % Update on each trial's result or based on batch group results
    
    dlTrainOptions('dlLearningRule') = 'BioDeltaRule'; % Delta rule with a basic change based on biophysical properties 
    dlTrainOptions('dlSimulationFlag') = 1; % If 0, will not run simulations (only for debugging purposes)
    dlTrainOptions('dlOutputLogFlag') = 0; % If 0, will not keep outputs
    dlTrainOptions('dlOfflineOutputGenerator') = 0; % If 1, will generate fake-random outputs (only for debugging purposes)
    
    dlTrainOptions('dlAdaptiveLambda') = 1; % Adaptive lambda parameter; recommended for long simulations.
    dlTrainOptions('dlLambdaCap') = 1.1; % Only if Adaptive lambda is active, recommended to set a upper-bound (UB) or ignore to use default UB (0.01).
    dlTrainOptions('dlExcludeDiverge') = 1; % Exclude non-optimals from model log
    dlTrainOptions('dlTrainExcludeList') = {'Stim', ['DeepISOM', ModelName, '->'], ['DeepIPV', ModelName, '->'], ['MidIPV', ModelName, '->'], ['supIPV', ModelName, '->'], ['supISOM', ModelName, '->']}; % Exclude populations from training
    
    argsPowSpectRatio1 = struct();
    argsPowSpectRatio2 = struct();
    argsPowSpectRatio3 = struct();
    argsPowSpectRatio4 = struct();
    
    argsPowSpectRatio1.lf1 = 2;
    argsPowSpectRatio1.hf1 = 6;
    argsPowSpectRatio2.lf1 = 8;
    argsPowSpectRatio2.hf1 = 14;
    
    argsPowSpectRatio3.lf1 = 15;
    argsPowSpectRatio3.hf1 = 30;
    argsPowSpectRatio4.lf1 = 40;
    argsPowSpectRatio4.hf1 = 90;
 
    dlTrainOptions('dlLambda') = 1e-9; % 1e-11(1) -> 1e-4 (4)
    dlTrainOptions('dlAdaptiveLambda') = 1; % Adaptive lambda parameter; recommended for long simulations.
    dlTrainOptions('dlUpdateMode') = 'trial';
    dlTrainOptions('dlLearningRule') = 'EnhancedDeltaRule';
    
    if model_size_id > 10
        dlTrainOptions('dlTrainExcludeList') = {'Stim', ['DeepISOM', ModelName, '->'], ['DeepIPV', ModelName, '->'], ['MidIPV', ModelName, '->'], ['supIPV', ModelName, '->'], ['supISOM', ModelName, '->']}; % Exclude PV+SOM
%     elseif model_size_id > 30
%         dlTrainOptions('dlTrainExcludeList') = {'Stim', ['deepIPV', ModelName, '->'], ['midIPV', ModelName, '->'], ['supIPV', ModelName, '->']}; % Exclude PV
%     elseif model_size_id > 20
%         dlTrainOptions('dlTrainExcludeList') = {'Stim', ['deepISOM', ModelName, '->'], ['supISOM', ModelName, '->']}; % Exclude SOM
%     elseif model_size_id > 10
%         dlTrainOptions('dlTrainExcludeList') = {'Stim', ['deepE', ModelName, '->'], ['midE', ModelName, '->'], ['supE', ModelName, '->']}; % Exclude Excitatories
    else
        dlTrainOptions('dlTrainExcludeList') = {'Stim'}; % Normal, all included.
    end
   
    dlTrainOptions('dlCheckpointLengthCap') = 14;
    dlTrainOptions('dlEpochs') = 1;
    dlTrainOptions('dlBatchs') = 100;
    
    dlTrainOptions('dlEnhancedMomentum') = 0.6;
    CheckCoeff = 2.01;

    if strcmpi(ResetOptimalError, 'on')
        m.dlResetTraining();
    end
    
    dlTrainOptions('dlCustomLog') = ["dlEPowerSpectrum", "dlEPowerSpectrum", "dlEPowerSpectrum", "dlEPowerSpectrum", "dlLFPaxLog", "dlAccuracyBastos2020Task"]; % Name of a function which is in the path
    dlTrainOptions('dlCustomLogArgs') = [argsPowSpectRatio1, argsPowSpectRatio2, argsPowSpectRatio3, argsPowSpectRatio4, argsPowSpectRatio1, argsPowSpectRatio1]; % Arguments of your custom function

    if tune_flag

        disp("--> Tune mode")
        for cnt = 1:1
        
            disp("----------U0----------");  
            m.dlOptimalError = 1e9;dlTrainOptions('dlExcludeDiverge') = 0;
            dlTrainOptions('dlCheckpointCoefficient') = CheckCoeff*1.2;
            m.dlTrain(TBdata.TrB, dlOutputParameters, TBdata.TrT, dlTrainOptions);
        
            disp("----------U1----------");
            m.dlErrorsLog = [m.dlErrorsLog, -1];  
            m.dlOptimalError = 1e9;dlTrainOptions('dlExcludeDiverge') = 0;
            dlTrainOptions('dlCheckpointCoefficient') = CheckCoeff*1.2;
            m.dlTrain(TBdata.TrB, dlOutputParameters, TBdata.TrT, dlTrainOptions);

        end

    else

        for cnt = 1:1
        
            disp("----------U0----------");  
            m.dlOptimalError = 1e9;dlTrainOptions('dlExcludeDiverge') = 0;
            dlTrainOptions('dlCheckpointCoefficient') = CheckCoeff*1.2;
            m.dlTrain(TBdata.TrB, dlOutputParameters, TBdata.TrT, dlTrainOptions);
            
            disp("----------A-----------");
            m.dlErrorsLog = [m.dlErrorsLog, -1];
            m.dlOptimalError = 1e9;dlTrainOptions('dlExcludeDiverge') = 1;
            dlTrainOptions('dlCheckpointCoefficient') = CheckCoeff;
            m.dlTrain(TBdata.B1, dlOutputParameters, TBdata.T1, dlTrainOptions);
            
            disp("----------U1----------");
            m.dlErrorsLog = [m.dlErrorsLog, -1];  
            m.dlOptimalError = 1e9;dlTrainOptions('dlExcludeDiverge') = 0;
            dlTrainOptions('dlCheckpointCoefficient') = CheckCoeff*1.2;
            m.dlTrain(TBdata.TrB, dlOutputParameters, TBdata.TrT, dlTrainOptions);
            
            disp("----------B-----------");
            m.dlErrorsLog = [m.dlErrorsLog, -1]; 
            m.dlOptimalError = 1e9;dlTrainOptions('dlExcludeDiverge') = 1;
            dlTrainOptions('dlCheckpointCoefficient') = CheckCoeff; 
            m.dlTrain(TBdata.B2, dlOutputParameters, TBdata.T2, dlTrainOptions);
            
            disp("----------U2----------");
            m.dlErrorsLog = [m.dlErrorsLog, -1]; 
            m.dlOptimalError = 1e9;dlTrainOptions('dlExcludeDiverge') = 0;
            dlTrainOptions('dlCheckpointCoefficient') = CheckCoeff*1.2; 
            m.dlTrain(TBdata.TrB, dlOutputParameters, TBdata.TrT, dlTrainOptions);
            
            disp("----------C----------");
            m.dlErrorsLog = [m.dlErrorsLog, -1]; 
            m.dlOptimalError = 1e9;dlTrainOptions('dlExcludeDiverge') = 1;
            dlTrainOptions('dlCheckpointCoefficient') = CheckCoeff;
            m.dlTrain(TBdata.B3, dlOutputParameters, TBdata.T3, dlTrainOptions);
            
            disp("----------U3----------");
            m.dlErrorsLog = [m.dlErrorsLog, -1]; 
            m.dlOptimalError = 1e9;dlTrainOptions('dlExcludeDiverge') = 0;
            dlTrainOptions('dlCheckpointCoefficient') = CheckCoeff*1.2;
            m.dlTrain(TBdata.TrB, dlOutputParameters, TBdata.TrT, dlTrainOptions);
        
        end

    end

    m.dlSave();

end