function [trialParams1, trialParams2, trialParams3] = dlDemoThreePattern(suffixLabel)

    trialParams1 = containers.Map();
    trialParams2 = containers.Map();
    trialParams3 = containers.Map();

    trialParams1('tspan') = [0 1000];
    trialParams2('tspan') = [0 1000];
    trialParams3('tspan') = [0 1000];

    trialParams1(['PreStimuliA', suffixLabel, 't1']) = 100;
    trialParams1(['PreStimuliA', suffixLabel, 't2']) = 250;    
    trialParams1(['PostStimuliA', suffixLabel, 't1']) = 300;
    trialParams1(['PostStimuliA', suffixLabel, 't2']) = 450;

    trialParams1(['PreStimuliB', suffixLabel, 't1']) = 250;
    trialParams1(['PreStimuliB', suffixLabel, 't2']) = 250;    
    trialParams1(['PostStimuliB', suffixLabel, 't1']) = 300;
    trialParams1(['PostStimuliB', suffixLabel, 't2']) = 450;
    
    trialParams1(['PreStimuliC', suffixLabel, 't1']) = 250;
    trialParams1(['PreStimuliC', suffixLabel, 't2']) = 250;    
    trialParams1(['PostStimuliC', suffixLabel, 't1']) = 300;
    trialParams1(['PostStimuliC', suffixLabel, 't2']) = 450;
    
    trialParams2(['PreStimuliA', suffixLabel, 't1']) = 250;
    trialParams2(['PreStimuliA', suffixLabel, 't2']) = 250;    
    trialParams2(['PostStimuliA', suffixLabel, 't1']) = 300;
    trialParams2(['PostStimuliA', suffixLabel, 't2']) = 450;

    trialParams2(['PreStimuliB', suffixLabel, 't1']) = 100;
    trialParams2(['PreStimuliB', suffixLabel, 't2']) = 250;    
    trialParams2(['PostStimuliB', suffixLabel, 't1']) = 300;
    trialParams2(['PostStimuliB', suffixLabel, 't2']) = 450;
    
    trialParams2(['PreStimuliC', suffixLabel, 't1']) = 250;
    trialParams2(['PreStimuliC', suffixLabel, 't2']) = 250;    
    trialParams2(['PostStimuliC', suffixLabel, 't1']) = 300;
    trialParams2(['PostStimuliC', suffixLabel, 't2']) = 450;
    
    trialParams3(['PreStimuliA', suffixLabel, 't1']) = 250;
    trialParams3(['PreStimuliA', suffixLabel, 't2']) = 250;    
    trialParams3(['PostStimuliA', suffixLabel, 't1']) = 300;
    trialParams3(['PostStimuliA', suffixLabel, 't2']) = 450;

    trialParams3(['PreStimuliB', suffixLabel, 't1']) = 250;
    trialParams3(['PreStimuliB', suffixLabel, 't2']) = 250;    
    trialParams3(['PostStimuliB', suffixLabel, 't1']) = 300;
    trialParams3(['PostStimuliB', suffixLabel, 't2']) = 450;
    
    trialParams3(['PreStimuliC', suffixLabel, 't1']) = 100;
    trialParams3(['PreStimuliC', suffixLabel, 't2']) = 250;    
    trialParams3(['PostStimuliC', suffixLabel, 't1']) = 300;
    trialParams3(['PostStimuliC', suffixLabel, 't2']) = 450;
    
end
