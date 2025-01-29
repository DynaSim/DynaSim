function trialParams = dlGradedInputs(duration, t0, t1, Iappmx, n)

    trialParams = cell(1);

    for i = 1:n
    
        trialParams{i} = containers.Map();

        trialParams{i}('tspan') = [0 duration];
        trialParams{i}('EXc_Iapp') = Iappmx*((i-1)/(n));
        trialParams{i}('INg_Iapp') = Iappmx*((i-1)/(n));
        trialParams{i}('INl_Iapp') = Iappmx*((i-1)/(n));


        trialParams{i}('EXc_omega') = 0.12;
        trialParams{i}('INg_omega') = 0.12;
        trialParams{i}('INl_omega') = 0.12;
    
        trialParams{i}('EXc_t0') = t0;
        trialParams{i}('EXc_t1') = t1;
        trialParams{i}('INg_t0') = t0;
        trialParams{i}('INg_t1') = t1;
        trialParams{i}('INl_t0') = t0;
        trialParams{i}('INl_t1') = t1;

    end
    
end
