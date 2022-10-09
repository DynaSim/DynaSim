function dlThreeCuesTaskPerformer(Currentsize, model_size_id)

    %%% Create model parameters struct
    ModelParametersPFC = struct();
    
    %%% Area PFC layer sizes (relative)
    ModelParametersPFC.NeSuperficial = ceil(0.3*Currentsize);
    ModelParametersPFC.NSomSuperficial = ceil(0.04*Currentsize);
    ModelParametersPFC.NPvSuperficial = ceil(0.07*Currentsize);
    ModelParametersPFC.NeMid = ceil(0.14*Currentsize);
    ModelParametersPFC.NSomMid = 0;
    ModelParametersPFC.NPvMid = ceil(0.04*Currentsize);
    ModelParametersPFC.NeDeep = ceil(0.45*Currentsize);
    ModelParametersPFC.NSomDeep = ceil(0.04*Currentsize);
    ModelParametersPFC.NPvDeep = ceil(0.01*Currentsize);
    
    ModelParametersPFC.Nin = 6;
    ModelParametersPFC.Nout = 6;
    ModelParametersPFC.NoiseRate = 4; % 6%
    ModelParametersPFC.Nstim = 3;
    
    %%% Call Laminar Cortex Constructor Functions
    dsCellPFC = dlLaminarCortexNetLWK(ModelParametersPFC, 'PFC'); % Laminar PFC model with specific parameters
    
    %%% Finalization
    dsModel = dsCellPFC;
    
    % Create DynaLearn Class 
    % Try to use this section only first time or If you have lost your file and
    % you want a new model.
    
    m = DynaLearn(dsModel, char("models/dlPredictiveCorticalCircuitModelLWK" + string(model_size_id)), 'mex'); % ~10 min or less, MEXGEN or < 20 sec, RAWGEN.
    m.dlSave(); % < 1sec
    
    % Load DynaLearn Class
    
    model_size_id = 1;
    m = DynaLearn(); % ~ 1sec
    m = m.dlLoad(char("models/dlPredictiveCorticalCircuitModelLWK" + string(model_size_id))); % ~ 10sec, New larger model; keeping track of its activity in Gamma/Beta **

    [trialParams1, trialParams2, trialParams3] = dlDemoThreePattern('xPFC');

    outputParams = [{'deepExPFC_V', 1:floor(ModelParametersPFC.NeDeep/3), [400 750] ...
        , 'afr'}; {'deepExPFC_V',ceil(ModelParametersPFC.NeDeep/3):floor(2*ModelParametersPFC.NeDeep/3), ...
        [400 750], 'afr'}; {'deepExPFC_V', ceil(2*ModelParametersPFC.NeDeep/3):ModelParametersPFC.NeDeep, [400 750], 'afr'}; ...
        {'supExPFC_V', 1:ModelParametersPFC.NeSuperficial, [300 800], 'afr'}; ...
        {'midExPFC_V', 1:ModelParametersPFC.NeMid, [300 800], 'afr'}; ...
        {'deepExPFC_V', 1:ModelParametersPFC.NeDeep, [300 800], 'afr'}; ...
        {'supIPVxPFC_V', 1:ModelParametersPFC.NPvSuperficial, [300 800], 'afr'}; ...
        {'midIPVxPFC_V', 1:ModelParametersPFC.NPvMid, [300 800], 'afr'}; ...
        {'deepIPVxPFC_V', 1:ModelParametersPFC.NPvDeep, [300 800], 'afr'}];
    
    targetParams1 = [{'BGPenalty', 4:7, 40, 0.4}; {'TotalSpikesPenalty', 4:7, 40, 0.2}; {'Compare', [1, 2], 0, 0.2}; {'Compare', [1, 3], 0, 0.2}]; % A 
    targetParams2 = [{'BGPenalty', 4:7, 40, 0.4}; {'TotalSpikesPenalty', 4:7, 40, 0.2}; {'Compare', [2, 1], 0, 0.2}; {'Compare', [2, 3], 0, 0.2}]; % B
    targetParams3 = [{'BGPenalty', 4:7, 40, 0.4}; {'TotalSpikesPenalty', 4:7, 40, 0.2}; {'Compare', [3, 1], 0, 0.2}; {'Compare', [3, 2], 0, 0.2}]; % C
    
    dlInputParameters = {trialParams1, trialParams2, trialParams3};
    dlTargetParameters = {targetParams1, targetParams2, targetParams3};
    dlOutputParameters = outputParams;
    
    TBdata = dlTrialBlockGenerator(dlInputParameters, dlTargetParameters, 25, 25);
    
    dlTrainOptions = containers.Map(); % Train options; MUST be a map data structure
    dlTrainOptions('dlEpochs') = 100; % % Number of epochs (A.K.A total iterations)
    dlTrainOptions('dlBatchs') = 3; % If a scenario requires the training to be based on a group parameter (e.g mean of errors) use a dlBatch > 1 and set update mode later to batch. 
    dlTrainOptions('dlLambda') = 1e-5; % Higher lambda means more changes based on error, lower may cause model to learn slower or nothing.
        
    dlTrainOptions('dlCheckpoint') = 'true'; % If current step's error is higher based on a threshold, reload last optimal state and continue from that point
    dlTrainOptions('dlCheckpointCoefficient') = 2.047; % A.K.A exploration rate
    dlTrainOptions('dlCheckpointLengthCap') = 7; % If more than 7 steps with no progress passed, return to last checkpoint.
    dlTrainOptions('dlUpdateMode') = 'batch'; % Update on each trial's result or based on batch group results
    
    dlTrainOptions('dlLearningRule') = 'BioDeltaRule'; % Delta rule with a basic change based on biophysical properties 
    dlTrainOptions('dlSimulationFlag') = 1; % If 0, will not run simulations (only for debugging purposes)
    dlTrainOptions('dlOutputLogFlag') = 1; % If 0, will not keep outputs
    dlTrainOptions('dlOfflineOutputGenerator') = 0; % If 1, will generate fake-random outputs (only for debugging purposes)
    
    dlTrainOptions('dlAdaptiveLambda') = 1; % Adaptive lambda parameter; recommended for long simulations.
    dlTrainOptions('dlLambdaCap') = 1.1; % Only if Adaptive lambda is active, recommended to set a upper-bound (UB) or ignore to use default UB (0.01).
    dlTrainOptions('dlExcludeDiverge') = 1; % Exclude non-optimals from model log
    dlTrainOptions('dlTrainExcludeList') = {'Stim'}; % Exclude populations from training
    
    dlTrainOptions('dlLambda') = 7e-4;
    dlTrainOptions('dlEpochs') = 10;
    dlTrainOptions('dlBatchs') = 3;
    
    argsPowSpectRatio = struct();
    argsNull = [];
    
    argsPowSpectRatio.lf1 = 7;
    argsPowSpectRatio.hf1 = 28;
    argsPowSpectRatio.lf2 = 42;
    argsPowSpectRatio.hf2 = 84;
    
    dlTrainOptions('dlCustomLog') = ["dlPowerSpectrumRatio", "dlAccuracyBastos2020Task"]; % Name of a function which is in the path
    dlTrainOptions('dlCustomLogArgs') = [argsPowSpectRatio, argsNull]; % Arguments of your custom function
    
    dlTrainOptions('dlLambda') = 1e-7; % 1e-11(1) -> 1e-4 (4)
    dlTrainOptions('dlAdaptiveLambda') = 0; % Adaptive lambda parameter; recommended for long simulations.
    dlTrainOptions('dlUpdateMode') = 'trial';
    dlTrainOptions('dlLearningRule') = 'BioDeltaRule';
    
    dlTrainOptions('dlTrainExcludeList') = {'Stimuli'};
    dlTrainOptions('dlCheckpointLengthCap') = 15;
    dlTrainOptions('dlEpochs') = 5;
    dlTrainOptions('dlBatchs') = 10;
    
    dlTrainOptions('dlEnhancedMomentum') = 0.7;
    CheckCoeff = 1.5;
    m.dlResetTraining();
    argsPSR = struct();
    
    argsPSR.lf1 = 12;
    argsPSR.hf1 = 30;
    argsPSR.lf2 = 40;
    argsPSR.hf2 = 100;
    
    dlTrainOptions('dlCustomLog') = ["dlPowerSpectrumRatio", "dlAccuracyBastos2020Task"]; % Name of a function which is in the path
    dlTrainOptions('dlCustomLogArgs') = [argsPSR, argsNull]; % Arguments of your custom function
    
    for cnt = 1:1
    
        disp("----------A-----------");
        dlTrainOptions('dlExcludeDiverge') = 1;
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

    m.dlSave();

end