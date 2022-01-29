%{
TODO: Update text
Purpose: Incorporate spike-timing-dependent plasticity (STDP learning) into any of your DynaSim spiking models.
"STDP_Isyn" can be incorporated into a Dynasim model as a connection mechanism between any two populations of neurons. 
This tutorial demonstrates how to use it in a simple toy model with two HH-type neurons (E1 and E2). 
The last plot that the script generates shows the synaptic weight increasing after each spike from E1 to E2.

The connection mechanism consists of AMPA and NMDA synapses with the "online implementation" of an STDP model with soft bounds described here: 
http://www.scholarpedia.org/article/Spike-timing_dependent_plasticity#Online_implementation_of_STDP_models
%}
%function [rate, raster, data] = spn_network(condition)

    if nargin < 1
        condition = 'spontaneous';
        % condition = 'DC';
        % condition = 'reduced_DC';
        % condition = 'AC';
        % condition = 'reduced_AC';
    end

    % Simulation time
    tOn = 0;
    tOff = 1500;
    tspan = [tOn,tOff];
    transient = 500;
    interval = diff(tspan);

    % Define equations of cell model
    eqns = {
    'dV/dt = @current'
    'V(0) = -65+5*randn(1,Npop)' % mV, Initial conditions
    'monitor @current'
    % 'threshold = 0' % mV, Spike threshold
    % 'monitor V.spikes(threshold)'
    % 'monitor iGABAa.functions'
    };
    % TIP: monitor all functions of a mechanism using: monitor MECHANISM.functions

    
    % Create DynaSim specification structure
    s = [];

    % Model parameters
    numPseudoColumns = 2;                       % 1; 2;
    nSPNs = numPseudoColumns*75;                % 1; numPseudoColumns*15; numPseudoColumns*1; numPseudoColumns*75 % Number of SPN cells in the pool
    multiplicity = [nSPNs,nSPNs];

    g_m = 1.3;

    g_l_D1 = 0.096;
    g_l_D2 = 0.1;

    g_cat_D1 = 0.018;
    g_cat_D2 = 0.025;

    % Connection probabilities
    conn_prob_gaba_D1toD1 = 0.26;
    conn_prob_gaba_D1toD2 = 0.06;
    conn_prob_gaba_D2toD1 = 0.27;
    conn_prob_gaba_D2toD2 = 0.36;
    conn_prob_gaba_M = [conn_prob_gaba_D1toD1,conn_prob_gaba_D1toD2;conn_prob_gaba_D2toD1,conn_prob_gaba_D2toD2];

    % Calculate Gaba conductance for each connection
    g_gaba_D1_lowDA = 0.5/multiplicity(1);
    g_gaba_D2_lowDA = 1.65*g_gaba_D1_lowDA*multiplicity(1)/multiplicity(2);
    f_gaba_D1_highDA = 1.3;
    f_gaba_D2_highDA = 0.5;

    g_gaba_D1 = g_gaba_D1_lowDA*f_gaba_D1_highDA;
    g_gaba_D2 = g_gaba_D2_lowDA*f_gaba_D2_highDA;
    cross_inhibition_factor = 1.5; % 1;

    % Calculate Tau Gaba for each connection
    median_tau_gaba = 30.4;      % ms
    sig_tau_gaba = 8.2;          % ms
    tau_gaba_D1 = median_tau_gaba + sig_tau_gaba*randn(1,multiplicity(1));
    tau_gaba_D2 = median_tau_gaba + sig_tau_gaba*randn(1,multiplicity(2));

    % Tau STD
    tauD_D1_highDA = 1030;    % from 275 % ms
    tauD_D2_highDA = 210;     % from 415

    % background Async mean input
    g_extInp = 0.00035;    % in mS/cm^2
    baseline_extInp = 30e3; % in spks/s

    % pfc input to SPN (excitatory beta input)

    g_pfcInp = g_extInp;
    DC_ref = 44e3; % in spks/s
    AC_ref = 8e3;

    % freq_ref = 25; % High beta
    % freq_ref = 20; % Mid beta
    freq_ref = 18; % Low beta
    % freq_ref = 10; % Alpha

    tOn_pfcInp = 500; % ms
    tOff_pfcInp = 1500; % ms
    latency_pfcInp = 40; % ms
    dutyPercentage_pfcInp = 50;

    switch condition
        case 'spontaneous'
            DC_pfcInp = 0;
            AC_pfcInp = 0;
            freq_pfcInp = 0;
            thresholddutycycle_pfcInp = -1;
        case 'DC'
            DC_pfcInp = DC_ref;
            AC_pfcInp = 0;
            freq_pfcInp = 0;
            thresholddutycycle_pfcInp = -1;
        case 'AC'
            DC_pfcInp = DC_ref;
            AC_pfcInp = AC_ref;
            freq_pfcInp = freq_ref;
            thresholddutycycle_pfcInp = cos(pi*dutyPercentage_pfcInp/100);
        case 'reduced_DC'
            DC_pfcInp = DC_ref/5;
            AC_pfcInp = 0;
            freq_pfcInp = 0;
            thresholddutycycle_pfcInp = -1;
        case 'reduced_AC'
            DC_pfcInp = DC_ref/5;
            AC_pfcInp = AC_ref/5;
            freq_pfcInp = freq_ref;
            thresholddutycycle_pfcInp = cos(pi*dutyPercentage_pfcInp/100);
        otherwise
            DC_pfcInp = 0;
            AC_pfcInp = 0;
            freq_pfcInp = 0;
            thresholddutycycle_pfcInp = -1;
    end

    % injected current for basal excitation
    %injectedCurrent_D1 = 0; % -2.6; % 0;
    %injectedCurrent_D2 = 0;

    % injected current for fIcurve analysis
    %MaxInjCurrent = 2.7;          % microA/cm^2
    %MinInjCurrent = -0.3;         % microA/cm^2
    %tOn_fIcurve = transient;
    %tOff_fIcurve = tOff;           % in ms
    %npoints_fIcurve = 11

    % Main structure
    s.populations(1).name = 'D1_SPN';
    s.populations(1).size = multiplicity(1);
    s.populations(1).equations = eqns;
    s.populations(1).mechanism_list = {'spn_iNa','spn_iK','spn_iLeak','spn_iM','spn_iCa','spn_CaBuffer','spn_iKca','spn_iPoisson','spn_ipfcPoisson'};
    s.populations(1).parameters = {'cm',1,'g_m',g_m,'g_cat',g_cat_D1,'g_poisson',g_extInp,'baseline_poisson',baseline_extInp,'g_pfc_poisson',g_pfcInp,'onset_pfc_poisson',tOn_pfcInp,'offset_pfc_poisson',tOff_pfcInp,'f_pfc_poisson',freq_pfcInp,'thresholddutycycle_pfc_poisson',thresholddutycycle_pfcInp,'latency_pfc_poisson',latency_pfcInp,'DC_pfc_poisson',DC_pfcInp,'AC_pfc_poisson',AC_pfcInp};
    s.populations(2).name = 'D2_SPN';
    s.populations(2).size = multiplicity(2);
    s.populations(2).equations = eqns;
    s.populations(2).mechanism_list = {'spn_iNa','spn_iK','spn_iLeak','spn_iM','spn_iCa','spn_CaBuffer','spn_iKca','spn_iPoisson','spn_ipfcPoisson'};
    s.populations(2).parameters = {'cm',1,'g_m',g_m,'g_cat',g_cat_D2,'g_poisson',g_extInp,'baseline_poisson',baseline_extInp,'g_pfc_poisson',g_pfcInp,'onset_pfc_poisson',tOn_pfcInp,'offset_pfc_poisson',tOff_pfcInp,'f_pfc_poisson',freq_pfcInp,'thresholddutycycle_pfc_poisson',thresholddutycycle_pfcInp,'latency_pfc_poisson',latency_pfcInp,'DC_pfc_poisson',DC_pfcInp,'AC_pfc_poisson',AC_pfcInp};

    s.connections(1).direction = 'D1_SPN->D1_SPN';
    s.connections(1).mechanism_list = {'spn_iGABA'};
    s.connections(1).parameters = {'conn_prob_gaba',conn_prob_gaba_D1toD1,'g_gaba',g_gaba_D1,'tau_gaba',tau_gaba_D1,'tauD_gaba',tauD_D1_highDA};
    s.connections(2).direction = 'D1_SPN->D2_SPN';
    s.connections(2).mechanism_list = {'spn_iGABA'};
    s.connections(2).parameters = {'conn_prob_gaba',conn_prob_gaba_D1toD2,'g_gaba',cross_inhibition_factor*g_gaba_D1,'tau_gaba', tau_gaba_D1,'tauD_gaba',tauD_D1_highDA};
    s.connections(3).direction = 'D2_SPN->D1_SPN';
    s.connections(3).mechanism_list = {'spn_iGABA'};
    s.connections(3).parameters = {'conn_prob_gaba',conn_prob_gaba_D2toD1,'g_gaba',cross_inhibition_factor*g_gaba_D2,'tau_gaba', tau_gaba_D2,'tauD_gaba',tauD_D2_highDA};
    s.connections(4).direction = 'D2_SPN->D2_SPN';
    s.connections(4).mechanism_list = {'spn_iGABA'};
    s.connections(4).parameters = {'conn_prob_gaba',conn_prob_gaba_D2toD2,'g_gaba',g_gaba_D2,'tau_gaba', tau_gaba_D2,'tauD_gaba',tauD_D2_highDA};

    % Simulate
    data = dsSimulate(s,'time_limits',tspan,'dt',0.05,'solver','rk4'); % ,'mex_flag',1);

    % Salva's raster plot
    raster{1} = computeRaster(data.time,data.D1_SPN_V);
    raster{2} = computeRaster(data.time,data.D2_SPN_V);
    tl = [tOn+transient, tOff]; % [tOn, tOff];
    rate = plotRaster(multiplicity,tl,raster);
%end
