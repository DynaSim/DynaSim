function [trialParams1, trialParams2, trialParams3] = dlDemoThreePattern()

    trialParams1 = containers.Map();
    trialParams2 = containers.Map();
    trialParams3 = containers.Map();

    trialParams1('tspan') = [0 1000];
    trialParams2('tspan') = [0 1000];
    trialParams3('tspan') = [0 1000];

    trialParams1('PreStimuliAxV4_t1') = 100;
    trialParams1('PreStimuliAxV4_t2') = 250;    
    trialParams1('PostStimuliAxV4_t1') = 300;
    trialParams1('PostStimuliAxV4_t2') = 450;

    trialParams1('PreStimuliBxV4_t1') = 250;
    trialParams1('PreStimuliBxV4_t2') = 250;    
    trialParams1('PostStimuliBxV4_t1') = 300;
    trialParams1('PostStimuliBxV4_t2') = 450;
    
    trialParams1('PreStimuliCxV4_t1') = 250;
    trialParams1('PreStimuliCxV4_t2') = 250;    
    trialParams1('PostStimuliCxV4_t1') = 300;
    trialParams1('PostStimuliCxV4_t2') = 450;
    
    trialParams2('PreStimuliAxV4_t1') = 250;
    trialParams2('PreStimuliAxV4_t2') = 250;    
    trialParams2('PostStimuliAxV4_t1') = 300;
    trialParams2('PostStimuliAxV4_t2') = 450;

    trialParams2('PreStimuliBxV4_t1') = 100;
    trialParams2('PreStimuliBxV4_t2') = 250;    
    trialParams2('PostStimuliBxV4_t1') = 300;
    trialParams2('PostStimuliBxV4_t2') = 450;
    
    trialParams2('PreStimuliCxV4_t1') = 250;
    trialParams2('PreStimuliCxV4_t2') = 250;    
    trialParams2('PostStimuliCxV4_t1') = 300;
    trialParams2('PostStimuliCxV4_t2') = 450;
    
    trialParams3('PreStimuliAxV4_t1') = 250;
    trialParams3('PreStimuliAxV4_t2') = 250;    
    trialParams3('PostStimuliAxV4_t1') = 300;
    trialParams3('PostStimuliAxV4_t2') = 450;

    trialParams3('PreStimuliBxV4_t1') = 250;
    trialParams3('PreStimuliBxV4_t2') = 250;    
    trialParams3('PostStimuliBxV4_t1') = 300;
    trialParams3('PostStimuliBxV4_t2') = 450;
    
    trialParams3('PreStimuliCxV4_t1') = 100;
    trialParams3('PreStimuliCxV4_t2') = 250;    
    trialParams3('PostStimuliCxV4_t1') = 300;
    trialParams3('PostStimuliCxV4_t2') = 450;
    
end
