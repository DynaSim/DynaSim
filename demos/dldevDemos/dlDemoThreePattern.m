function [trialParams1, trialParams2, trialParams3] = dlDemoThreePattern()

    trialParams1 = containers.Map();
    trialParams2 = containers.Map();
    trialParams3 = containers.Map();

    trialParams1('tspan') = [0 1000];
    trialParams2('tspan') = [0 1000];
    trialParams3('tspan') = [0 1000];

    trialParams1('EAxV4_t1') = 100;
    trialParams1('EAxV4_t2') = 250;    
    trialParams1('IAxV4_t1') = 300;
    trialParams1('IAxV4_t2') = 450;

    trialParams1('EBxV4_t1') = 250;
    trialParams1('EBxV4_t2') = 250;    
    trialParams1('IBxV4_t1') = 300;
    trialParams1('IBxV4_t2') = 450;
    
    trialParams1('ECxV4_t1') = 250;
    trialParams1('ECxV4_t2') = 250;    
    trialParams1('ICxV4_t1') = 300;
    trialParams1('ICxV4_t2') = 450;
    
    trialParams2('EAxV4_t1') = 250;
    trialParams2('EAxV4_t2') = 250;    
    trialParams2('IAxV4_t1') = 300;
    trialParams2('IAxV4_t2') = 450;

    trialParams2('EBxV4_t1') = 100;
    trialParams2('EBxV4_t2') = 250;    
    trialParams2('IBxV4_t1') = 300;
    trialParams2('IBxV4_t2') = 450;
    
    trialParams2('ECxV4_t1') = 250;
    trialParams2('ECxV4_t2') = 250;    
    trialParams2('ICxV4_t1') = 300;
    trialParams2('ICxV4_t2') = 450;
    
    trialParams3('EAxV4_t1') = 250;
    trialParams3('EAxV4_t2') = 250;    
    trialParams3('IAxV4_t1') = 300;
    trialParams3('IAxV4_t2') = 450;

    trialParams3('EBxV4_t1') = 250;
    trialParams3('EBxV4_t2') = 250;    
    trialParams3('IBxV4_t1') = 300;
    trialParams3('IBxV4_t2') = 450;
    
    trialParams3('ECxV4_t1') = 100;
    trialParams3('ECxV4_t2') = 250;    
    trialParams3('ICxV4_t1') = 300;
    trialParams3('ICxV4_t2') = 450;
    
end
