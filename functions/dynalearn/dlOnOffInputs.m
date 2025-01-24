function [trialParams1, trialParams2] = dlOnOffInputs(duration, t0, t1, Iapp)

    trialParams1 = containers.Map();
    trialParams2 = containers.Map();

    trialParams1('tspan') = [0 duration];
    trialParams2('tspan') = [0 duration];

    % trialParams1('EXc_omega') = 0;
    % trialParams2('EXc_omega') = 0;

    trialParams1('EXc_Iapp') = 0;
    trialParams2('EXc_Iapp') = Iapp;
    
    trialParams1('EXc_t0') = t0;
    trialParams2('EXc_t0') = t0;

    trialParams1('EXc_t1') = t1;
    trialParams2('EXc_t1') = t1;
    
end
