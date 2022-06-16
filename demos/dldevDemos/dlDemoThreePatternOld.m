function [trialParams1, trialParams2, trialParams3] = dlDemoThreePatternOld()

    trialParams1 = containers.Map();
    trialParams2 = containers.Map();
    trialParams3 = containers.Map();

    trialParams1('tspan') = [0 500];
    trialParams2('tspan') = [0 500];
    trialParams3('tspan') = [0 500];

    dc_poisson = 7e5;

    trialParams1('SA1_ctx_iPoisson_DC_poisson') = dc_poisson;
    trialParams1('SA2_ctx_iPoisson_DC_poisson') = dc_poisson;
    trialParams1('SB1_ctx_iPoisson_DC_poisson') = 0;
    trialParams1('SB2_ctx_iPoisson_DC_poisson') = 0;
    trialParams1('SC1_ctx_iPoisson_DC_poisson') = 0;
    trialParams1('SC2_ctx_iPoisson_DC_poisson') = 0;

    trialParams1('SA1_ctx_iPoisson_onset_poisson') = 150;
    trialParams1('SA1_ctx_iPoisson_offset_poisson') = 250;
    trialParams1('SA2_ctx_iPoisson_onset_poisson') = 250;
    trialParams1('SA2_ctx_iPoisson_offset_poisson') = 350;

    trialParams1('SB1_ctx_iPoisson_onset_poisson') = 0;
    trialParams1('SB1_ctx_iPoisson_offset_poisson') = 0;
    trialParams1('SB2_ctx_iPoisson_onset_poisson') = 0;
    trialParams1('SB2_ctx_iPoisson_offset_poisson') = 0;

    trialParams1('SC1_ctx_iPoisson_onset_poisson') = 0;
    trialParams1('SC1_ctx_iPoisson_offset_poisson') = 0;
    trialParams1('SC2_ctx_iPoisson_onset_poisson') = 0;
    trialParams1('SC2_ctx_iPoisson_offset_poisson') = 0;

    trialParams2('SA1_ctx_iPoisson_DC_poisson') = dc_poisson;
    trialParams2('SA2_ctx_iPoisson_DC_poisson') = dc_poisson;
    trialParams2('SB1_ctx_iPoisson_DC_poisson') = dc_poisson;
    trialParams2('SB2_ctx_iPoisson_DC_poisson') = dc_poisson;
    trialParams2('SC1_ctx_iPoisson_DC_poisson') = dc_poisson;
    trialParams2('SC2_ctx_iPoisson_DC_poisson') = dc_poisson;

    trialParams2('SA1_ctx_iPoisson_onset_poisson') = 250;
    trialParams2('SA1_ctx_iPoisson_offset_poisson') = 250;
    trialParams2('SA2_ctx_iPoisson_onset_poisson') = 350;
    trialParams2('SA2_ctx_iPoisson_offset_poisson') = 350;

    trialParams2('SB1_ctx_iPoisson_onset_poisson') = 150;
    trialParams2('SB1_ctx_iPoisson_offset_poisson') = 250;
    trialParams2('SB2_ctx_iPoisson_onset_poisson') = 250;
    trialParams2('SB2_ctx_iPoisson_offset_poisson') = 350;

    trialParams2('SC1_ctx_iPoisson_onset_poisson') = 250;
    trialParams2('SC1_ctx_iPoisson_offset_poisson') = 250;
    trialParams2('SC2_ctx_iPoisson_onset_poisson') = 350;
    trialParams2('SC2_ctx_iPoisson_offset_poisson') = 350;

    trialParams3('SA1_ctx_iPoisson_DC_poisson') = dc_poisson;
    trialParams3('SA2_ctx_iPoisson_DC_poisson') = dc_poisson;
    trialParams3('SB1_ctx_iPoisson_DC_poisson') = dc_poisson;
    trialParams3('SB2_ctx_iPoisson_DC_poisson') = dc_poisson;
    trialParams3('SC1_ctx_iPoisson_DC_poisson') = dc_poisson;
    trialParams3('SC2_ctx_iPoisson_DC_poisson') = dc_poisson;

    trialParams3('SA1_ctx_iPoisson_onset_poisson') = 250;
    trialParams3('SA1_ctx_iPoisson_offset_poisson') = 250;
    trialParams3('SA2_ctx_iPoisson_onset_poisson') = 350;
    trialParams3('SA2_ctx_iPoisson_offset_poisson') = 350;

    trialParams3('SB1_ctx_iPoisson_onset_poisson') = 250;
    trialParams3('SB1_ctx_iPoisson_offset_poisson') = 250;
    trialParams3('SB2_ctx_iPoisson_onset_poisson') = 350;
    trialParams3('SB2_ctx_iPoisson_offset_poisson') = 350;

    trialParams3('SC1_ctx_iPoisson_onset_poisson') = 150;
    trialParams3('SC1_ctx_iPoisson_offset_poisson') = 250;
    trialParams3('SC2_ctx_iPoisson_onset_poisson') = 250;
    trialParams3('SC2_ctx_iPoisson_offset_poisson') = 350;
    
    end
