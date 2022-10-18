function [trialParams1, trialParams2, trialParams3] = dlDemoThreePattern(suffixLabel)

    trialParams1 = containers.Map();
    trialParams2 = containers.Map();
    trialParams3 = containers.Map();

    trialParams1('tspan') = [0 750];
    trialParams2('tspan') = [0 750];
    trialParams3('tspan') = [0 750];

    trialParams1(['PreStimuliA', suffixLabel, '_t1']) = 100;
    trialParams1(['PreStimuliA', suffixLabel, '_t2']) = 250;    
    trialParams1(['PostStimuliA', suffixLabel, '_t1']) = 300;
    trialParams1(['PostStimuliA', suffixLabel, '_t2']) = 450;

    trialParams1(['PreStimuliB', suffixLabel, '_t1']) = 250;
    trialParams1(['PreStimuliB', suffixLabel, '_t2']) = 250;    
    trialParams1(['PostStimuliB', suffixLabel, '_t1']) = 300;
    trialParams1(['PostStimuliB', suffixLabel, '_t2']) = 450;
    
    trialParams1(['PreStimuliC', suffixLabel, '_t1']) = 250;
    trialParams1(['PreStimuliC', suffixLabel, '_t2']) = 250;    
    trialParams1(['PostStimuliC', suffixLabel, '_t1']) = 300;
    trialParams1(['PostStimuliC', suffixLabel, '_t2']) = 450;
    
    trialParams2(['PreStimuliA', suffixLabel, '_t1']) = 250;
    trialParams2(['PreStimuliA', suffixLabel, '_t2']) = 250;    
    trialParams2(['PostStimuliA', suffixLabel, '_t1']) = 300;
    trialParams2(['PostStimuliA', suffixLabel, '_t2']) = 450;

    trialParams2(['PreStimuliB', suffixLabel, '_t1']) = 100;
    trialParams2(['PreStimuliB', suffixLabel, '_t2']) = 250;    
    trialParams2(['PostStimuliB', suffixLabel, '_t1']) = 300;
    trialParams2(['PostStimuliB', suffixLabel, '_t2']) = 450;
    
    trialParams2(['PreStimuliC', suffixLabel, '_t1']) = 250;
    trialParams2(['PreStimuliC', suffixLabel, '_t2']) = 250;    
    trialParams2(['PostStimuliC', suffixLabel, '_t1']) = 300;
    trialParams2(['PostStimuliC', suffixLabel, '_t2']) = 450;
    
    trialParams3(['PreStimuliA', suffixLabel, '_t1']) = 250;
    trialParams3(['PreStimuliA', suffixLabel, '_t2']) = 250;    
    trialParams3(['PostStimuliA', suffixLabel, '_t1']) = 300;
    trialParams3(['PostStimuliA', suffixLabel, '_t2']) = 450;

    trialParams3(['PreStimuliB', suffixLabel, '_t1']) = 250;
    trialParams3(['PreStimuliB', suffixLabel, '_t2']) = 250;    
    trialParams3(['PostStimuliB', suffixLabel, '_t1']) = 300;
    trialParams3(['PostStimuliB', suffixLabel, '_t2']) = 450;
    
    trialParams3(['PreStimuliC', suffixLabel, '_t1']) = 100;
    trialParams3(['PreStimuliC', suffixLabel, '_t2']) = 250;    
    trialParams3(['PostStimuliC', suffixLabel, '_t1']) = 300;
    trialParams3(['PostStimuliC', suffixLabel, '_t2']) = 450;
    
end
