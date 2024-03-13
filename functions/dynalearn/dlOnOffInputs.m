function [trialParams1, trialParams2] = dlOnOffInputs(duration, t0, t1, Iapp)

    trialParams1 = containers.Map();
    trialParams2 = containers.Map();

    trialParams1('tspan') = [0 duration];
    trialParams2('tspan') = [0 duration];

    trialParams1('ES_Iapp') = 0;
    trialParams2('ES_Iapp') = Iapp;
    
    trialParams1('ES_t0') = t0;
    trialParams2('ES_t0') = t0;

    trialParams1('ES_t1') = t1;
    trialParams2('ES_t1') = t1;
    
end
