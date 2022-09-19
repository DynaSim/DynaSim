function [trialParams1, trialParams2, trialParams3] = dlDemoThreePattern()

    trialParams1 = containers.Map();
    trialParams2 = containers.Map();
    trialParams3 = containers.Map();

    trialParams1('tspan') = [0 550];
    trialParams2('tspan') = [0 550];
    trialParams3('tspan') = [0 550];

    trialParams1('EAxV4_t1') = 70;
    trialParams1('EAxV4_t2') = 210;    
    trialParams1('IAxV4_t1') = 280;
    trialParams1('IAxV4_t2') = 430;

    trialParams1('EBxV4_t1') = 210;
    trialParams1('EBxV4_t2') = 210;    
    trialParams1('IBxV4_t1') = 280;
    trialParams1('IBxV4_t2') = 430;
    
    trialParams1('ECxV4_t1') = 210;
    trialParams1('ECxV4_t2') = 210;    
    trialParams1('ICxV4_t1') = 280;
    trialParams1('ICxV4_t2') = 430;
    
    trialParams2('EAxV4_t1') = 210;
    trialParams2('EAxV4_t2') = 210;    
    trialParams2('IAxV4_t1') = 280;
    trialParams2('IAxV4_t2') = 430;

    trialParams2('EBxV4_t1') = 70;
    trialParams2('EBxV4_t2') = 210;    
    trialParams2('IBxV4_t1') = 280;
    trialParams2('IBxV4_t2') = 430;
    
    trialParams2('ECxV4_t1') = 210;
    trialParams2('ECxV4_t2') = 210;    
    trialParams2('ICxV4t1') = 280;
    trialParams2('ICxV4_t2') = 430;
    
    trialParams3('EAxV4_t1') = 210;
    trialParams3('EAxV4_t2') = 210;    
    trialParams3('IAxV4_t1') = 280;
    trialParams3('IAxV4_t2') = 430;

    trialParams3('EBxV4_t1') = 210;
    trialParams3('EBxV4_t2') = 210;    
    trialParams3('IBxV4_t1') = 280;
    trialParams3('IBxV4_t2') = 430;
    
    trialParams3('ECxV4_t1') = 70;
    trialParams3('ECxV4_t2') = 210;    
    trialParams3('ICxV4_t1') = 280;
    trialParams3('ICxV4_t2') = 430;
    
end
