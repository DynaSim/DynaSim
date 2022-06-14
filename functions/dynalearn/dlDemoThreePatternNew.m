function [trialParams1, trialParams2, trialParams3] = dlDemoThreePatternNew()

    trialParams1 = containers.Map();
    trialParams2 = containers.Map();
    trialParams3 = containers.Map();

    trialParams1('tspan') = [0 500];
    trialParams2('tspan') = [0 500];
    trialParams3('tspan') = [0 500];

    trialParams1('IO_SA1_t1') = 200;
    trialParams1('IO_SA1_t2') = 300;    
    trialParams1('IO_SA2_t1') = 300;
    trialParams1('IO_SA2_t2') = 400;

    trialParams1('IO_SB1_t1') = 200;
    trialParams1('IO_SB1_t2') = 300;    
    trialParams1('IO_SB2_t1') = 300;
    trialParams1('IO_SB2_t2') = 300;
    
    trialParams1('IO_SC1_t1') = 200;
    trialParams1('IO_SC1_t2') = 300;    
    trialParams1('IO_SC2_t1') = 300;
    trialParams1('IO_SC2_t2') = 300;
    
    trialParams2('IO_SA1_t1') = 200;
    trialParams2('IO_SA1_t2') = 300;    
    trialParams2('IO_SA2_t1') = 300;
    trialParams2('IO_SA2_t2') = 300;

    trialParams2('IO_SB1_t1') = 200;
    trialParams2('IO_SB1_t2') = 300;    
    trialParams2('IO_SB2_t1') = 300;
    trialParams2('IO_SB2_t2') = 400;
    
    trialParams2('IO_SC1_t1') = 200;
    trialParams2('IO_SC1_t2') = 300;    
    trialParams2('IO_SC2_t1') = 300;
    trialParams2('IO_SC2_t2') = 300;
    
    trialParams3('IO_SA1_t1') = 200;
    trialParams3('IO_SA1_t2') = 300;    
    trialParams3('IO_SA2_t1') = 300;
    trialParams3('IO_SA2_t2') = 300;

    trialParams3('IO_SB1_t1') = 200;
    trialParams3('IO_SB1_t2') = 300;    
    trialParams3('IO_SB2_t1') = 300;
    trialParams3('IO_SB2_t2') = 300;
    
    trialParams3('IO_SC1_t1') = 200;
    trialParams3('IO_SC1_t2') = 300;    
    trialParams3('IO_SC2_t1') = 300;
    trialParams3('IO_SC2_t2') = 400;
    
end
