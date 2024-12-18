function [trialParams1, trialParams2, trialParams3] = dlDemoThreePattern(suffixLabel, tspan)

    trialParams1 = containers.Map();
    trialParams2 = containers.Map();
    trialParams3 = containers.Map();

    trialParams1('tspan') = [0 tspan];
    trialParams2('tspan') = [0 tspan];
    trialParams3('tspan') = [0 tspan];

    trialParams1(['IOSA1', suffixLabel, '_t1']) = 100;
    trialParams1(['IOSA1', suffixLabel, '_t2']) = 250;    
    trialParams1(['IOSA2', suffixLabel, '_t1']) = 300;
    trialParams1(['IOSA2', suffixLabel, '_t2']) = 450;

    trialParams1(['IOSB1', suffixLabel, '_t1']) = 250;
    trialParams1(['IOSB1', suffixLabel, '_t2']) = 250;    
    trialParams1(['IOSB2', suffixLabel, '_t1']) = 300;
    trialParams1(['IOSB2', suffixLabel, '_t2']) = 450;
    
    trialParams1(['IOSC1', suffixLabel, '_t1']) = 250;
    trialParams1(['IOSC1', suffixLabel, '_t2']) = 250;    
    trialParams1(['IOSC2', suffixLabel, '_t1']) = 300;
    trialParams1(['IOSC2', suffixLabel, '_t2']) = 450;
    
    trialParams2(['IOSA1', suffixLabel, '_t1']) = 250;
    trialParams2(['IOSA1', suffixLabel, '_t2']) = 250;    
    trialParams2(['IOSA2', suffixLabel, '_t1']) = 300;
    trialParams2(['IOSA2', suffixLabel, '_t2']) = 450;

    trialParams2(['IOSB1', suffixLabel, '_t1']) = 100;
    trialParams2(['IOSB1', suffixLabel, '_t2']) = 250;    
    trialParams2(['IOSB2', suffixLabel, '_t1']) = 300;
    trialParams2(['IOSB2', suffixLabel, '_t2']) = 450;
    
    trialParams2(['IOSC1', suffixLabel, '_t1']) = 250;
    trialParams2(['IOSC1', suffixLabel, '_t2']) = 250;    
    trialParams2(['IOSC2', suffixLabel, '_t1']) = 300;
    trialParams2(['IOSC2', suffixLabel, '_t2']) = 450;
    
    trialParams3(['IOSA1', suffixLabel, '_t1']) = 250;
    trialParams3(['IOSA1', suffixLabel, '_t2']) = 250;    
    trialParams3(['IOSA2', suffixLabel, '_t1']) = 300;
    trialParams3(['IOSA2', suffixLabel, '_t2']) = 450;

    trialParams3(['IOSB1', suffixLabel, '_t1']) = 250;
    trialParams3(['IOSB1', suffixLabel, '_t2']) = 250;    
    trialParams3(['IOSB2', suffixLabel, '_t1']) = 300;
    trialParams3(['IOSB2', suffixLabel, '_t2']) = 450;
    
    trialParams3(['IOSC1', suffixLabel, '_t1']) = 100;
    trialParams3(['IOSC1', suffixLabel, '_t2']) = 250;    
    trialParams3(['IOSC2', suffixLabel, '_t1']) = 300;
    trialParams3(['IOSC2', suffixLabel, '_t2']) = 450;
    
end
