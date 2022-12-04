function m = dlThreeCuesTaskPerformer(Currentsize, model_size_idc)

    %%% Create model parameters struct
    ModelParametersPFC = struct();
    
    %%% Area PFC layer sizes (relative)
    ModelParametersPFC.NeSuperficial = ceil(0.3*Currentsize);
    ModelParametersPFC.NSomSuperficial = ceil(0.03*Currentsize);
    ModelParametersPFC.NPvSuperficial = ceil(0.03*Currentsize);
    ModelParametersPFC.NeMid = ceil(0.14*Currentsize);
    ModelParametersPFC.NSomMid = 0;
    ModelParametersPFC.NPvMid = ceil(0.03*Currentsize);
    ModelParametersPFC.NeDeep = ceil(0.45*Currentsize);
    ModelParametersPFC.NSomDeep = ceil(0.04*Currentsize);
    ModelParametersPFC.NPvDeep = ceil(0.01*Currentsize);
    
    ModelParametersPFC.Nin = 6;
    ModelParametersPFC.Nout = 6;
    ModelParametersPFC.NoiseRate = noise_rate; % 10%
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

    m = DynaLearn(); % ~ 1sec
    m = m.dlLoad(char("models/dlPredictiveCorticalCircuitModelLWK" + string(model_size_id))); % ~ 10sec, New larger model; keeping track of its activity in Gamma/Beta **

    [trialParams1, trialParams2, trialParams3] = dlDemoThreePattern('xPFC');
    
    outputParams = [{'deepExPFC_V', 1:floor(ModelParametersPFC.NeDeep/3), [200 400] ...
        , 'afr'}; {'deepExPFC_V',ceil(ModelParametersPFC.NeDeep/3):floor(2*ModelParametersPFC.NeDeep/3), ...
        [200 400], 'afr'}; {'deepExPFC_V', ceil(2*ModelParametersPFC.NeDeep/3):ModelParametersPFC.NeDeep, [200 400], 'afr'}; ...
        {'supExPFC_V', 1:ModelParametersPFC.NeSuperficial, [50 700], 'afr'}; ...
        {'midExPFC_V', 1:ModelParametersPFC.NeMid, [50 700], 'afr'}; ...
        {'deepExPFC_V', 1:ModelParametersPFC.NeDeep, [50 700], 'afr'}; ...
        {'supIPVxPFC_V', 1:ModelParametersPFC.NPvSuperficial, [50 700], 'afr'}; ...
        {'midIPVxPFC_V', 1:ModelParametersPFC.NPvMid, [50 700], 'afr'}; ...
        {'deepIPVxPFC_V', 1:ModelParametersPFC.NPvDeep, [50 700], 'afr'}];
    
    targetParams1 = [{'EPenalty', 4:9, 200, 0.01}; {'Compare', [1, 2], 0, 5.5}; {'Compare', [1, 3], 0, 5.5}]; % A 
    targetParams2 = [{'EPenalty', 4:9, 200, 0.01}; {'Compare', [2, 1], 0, 5.5}; {'Compare', [2, 3], 0, 5.5}]; % B
    targetParams3 = [{'EPenalty', 4:9, 200, 0.01}; {'Compare', [3, 1], 0, 5.5}; {'Compare', [3, 2], 0, 5.5}]; % C
    
    dlInputParameters = {trialParams1, trialParams2, trialParams3};
    dlTargetParameters = {targetParams1, targetParams2, targetParams3};
    dlOutputParameters = outputParams;
    
    TBdata = dlTrialBlockGenerator(dlInputParameters, dlTargetParameters, 100, 100);
    
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
    dlTrainOptions('dlTrainExcludeList') = {'Stim', 'deepISOMxPFC->', 'deepIPVxPFC->', 'midIPVxPFC->', 'supIPVxPFC->', 'supISOMxPFC->'}; % Exclude populations from training
    
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
    
    if model_size_id > 20
        dlTrainOptions('dlTrainExcludeList') = {'Stim', 'deepISOMxPFC->', 'deepIPVxPFC->', 'midIPVxPFC->', 'supIPVxPFC->', 'supISOMxPFC->'}; % Exclude PV+SOM
%     elseif model_size_id > 30
%         dlTrainOptions('dlTrainExcludeList') = {'Stim', 'deepIPVxPFC->', 'midIPVxPFC->', 'supIPVxPFC->'}; % Exclude PV
%     elseif model_size_id > 20
%         dlTrainOptions('dlTrainExcludeList') = {'Stim', 'deepISOMxPFC->', 'supISOMxPFC->'}; % Exclude SOM
    elseif model_size_id > 10
        dlTrainOptions('dlTrainExcludeList') = {'Stim', 'deepExPFC->', 'midExPFC->', 'supExPFC->'}; % Exclude Excitatories
    else
        dlTrainOptions('dlTrainExcludeList') = {'Stim'}; % Normal, all included.
    end
   
    dlTrainOptions('dlCheckpointLengthCap') = 14;
    dlTrainOptions('dlEpochs') = 1;
    dlTrainOptions('dlBatchs') = 100;
    
    dlTrainOptions('dlEnhancedMomentum') = 0.6;
    CheckCoeff = 2.01;
    m.dlResetTraining();
    
    dlTrainOptions('dlCustomLog') = ["dlEPowerSpectrum", "dlEPowerSpectrum", "dlEPowerSpectrum", "dlEPowerSpectrum", "dlLFPaxLog", "dlAccuracyBastos2020Task"]; % Name of a function which is in the path
    dlTrainOptions('dlCustomLogArgs') = [argsPowSpectRatio1, argsPowSpectRatio2, argsPowSpectRatio3, argsPowSpectRatio4, argsPowSpectRatio1, argsPowSpectRatio1]; % Arguments of your custom function

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

    m.dlSave();

end