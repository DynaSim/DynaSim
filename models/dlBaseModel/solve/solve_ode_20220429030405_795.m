function [T,E1_V,E1_iNa_m,E1_iNa_h,E1_iK_n,I1_V,I1_iNa_m,I1_iNa_h,I1_iK_n,S1_V,S1_iNa_m,S1_iNa_h,S1_iK_n,S2_V,S2_iNa_m,S2_iNa_h,S2_iK_n,I1_E1_iAMPActx_s,E1_E1_iAMPActx_s,E1_I1_iGABActx_s,I1_I1_iGABActx_s,S2_S1_iAMPActx_s,S1_S1_iAMPActx_s,S1_S2_iAMPActx_s,S2_S2_iAMPActx_s,E1_S1_iAMPActx_s,E1_S2_iAMPActx_s,E1_ctx_iPoisson_g_poisson,E1_ctx_iPoisson_I_poisson,I1_ctx_iPoisson_g_poisson,I1_ctx_iPoisson_I_poisson,S1_ctx_iPoisson_g_poisson,S1_ctx_iPoisson_I_poisson,S2_ctx_iPoisson_g_poisson,S2_ctx_iPoisson_I_poisson,E1_ctx_iPoisson_s_poisson,I1_ctx_iPoisson_s_poisson,S1_ctx_iPoisson_s_poisson,S2_ctx_iPoisson_s_poisson]=solve_ode

% ------------------------------------------------------------
% Parameters:
% ------------------------------------------------------------
params = load('params.mat','p');
p = params.p;
downsample_factor=p.downsample_factor;
dt=p.dt;
T=(p.tspan(1):dt:p.tspan(2))';
ntime=length(T);
nsamp=length(1:downsample_factor:ntime);

% seed the random number generator
rng(p.random_seed);

% ------------------------------------------------------------
% Fixed variables:
% ------------------------------------------------------------
E1_ctx_iPoisson_s_poisson = getPoissonDutyCycleGating(p.E1_ctx_iPoisson_onset_poisson,p.E1_ctx_iPoisson_offset_poisson,p.E1_ctx_iPoisson_latency_poisson,p.E1_ctx_iPoisson_f_poisson,p.E1_ctx_iPoisson_fsigma_poisson,p.E1_ctx_iPoisson_phi_poisson,p.E1_ctx_iPoisson_connsigma_poisson,p.E1_ctx_iPoisson_baseline_poisson,p.E1_ctx_iPoisson_DC_poisson,p.E1_ctx_iPoisson_AC_poisson,p.E1_ctx_iPoisson_slopedutycycle_poisson,p.E1_ctx_iPoisson_thresholddutycycle_poisson,p.E1_ctx_iPoisson_tau_poisson,p.E1_ctx_iPoisson_kick_poisson,T,p.E1_Npop);
I1_ctx_iPoisson_s_poisson = getPoissonDutyCycleGating(p.I1_ctx_iPoisson_onset_poisson,p.I1_ctx_iPoisson_offset_poisson,p.I1_ctx_iPoisson_latency_poisson,p.I1_ctx_iPoisson_f_poisson,p.I1_ctx_iPoisson_fsigma_poisson,p.I1_ctx_iPoisson_phi_poisson,p.I1_ctx_iPoisson_connsigma_poisson,p.I1_ctx_iPoisson_baseline_poisson,p.I1_ctx_iPoisson_DC_poisson,p.I1_ctx_iPoisson_AC_poisson,p.I1_ctx_iPoisson_slopedutycycle_poisson,p.I1_ctx_iPoisson_thresholddutycycle_poisson,p.I1_ctx_iPoisson_tau_poisson,p.I1_ctx_iPoisson_kick_poisson,T,p.I1_Npop);
S1_ctx_iPoisson_s_poisson = getPoissonDutyCycleGating(p.S1_ctx_iPoisson_onset_poisson,p.S1_ctx_iPoisson_offset_poisson,p.S1_ctx_iPoisson_latency_poisson,p.S1_ctx_iPoisson_f_poisson,p.S1_ctx_iPoisson_fsigma_poisson,p.S1_ctx_iPoisson_phi_poisson,p.S1_ctx_iPoisson_connsigma_poisson,p.S1_ctx_iPoisson_baseline_poisson,p.S1_ctx_iPoisson_DC_poisson,p.S1_ctx_iPoisson_AC_poisson,p.S1_ctx_iPoisson_slopedutycycle_poisson,p.S1_ctx_iPoisson_thresholddutycycle_poisson,p.S1_ctx_iPoisson_tau_poisson,p.S1_ctx_iPoisson_kick_poisson,T,p.S1_Npop);
S2_ctx_iPoisson_s_poisson = getPoissonDutyCycleGating(p.S2_ctx_iPoisson_onset_poisson,p.S2_ctx_iPoisson_offset_poisson,p.S2_ctx_iPoisson_latency_poisson,p.S2_ctx_iPoisson_f_poisson,p.S2_ctx_iPoisson_fsigma_poisson,p.S2_ctx_iPoisson_phi_poisson,p.S2_ctx_iPoisson_connsigma_poisson,p.S2_ctx_iPoisson_baseline_poisson,p.S2_ctx_iPoisson_DC_poisson,p.S2_ctx_iPoisson_AC_poisson,p.S2_ctx_iPoisson_slopedutycycle_poisson,p.S2_ctx_iPoisson_thresholddutycycle_poisson,p.S2_ctx_iPoisson_tau_poisson,p.S2_ctx_iPoisson_kick_poisson,T,p.S2_Npop);

% ------------------------------------------------------------
% Initial conditions:
% ------------------------------------------------------------
t=0; k=1;

% STATE_VARIABLES:
E1_V_last =  -rand(1, p.E1_Npop)*20;
E1_V = zeros(nsamp,p.E1_Npop);
E1_V(1,:) = E1_V_last;
E1_iNa_m_last = p.E1_iNa_m_IC+p.E1_iNa_IC_noise*rand(1,p.E1_Npop);
E1_iNa_m = zeros(nsamp,p.E1_Npop);
E1_iNa_m(1,:) = E1_iNa_m_last;
E1_iNa_h_last = p.E1_iNa_h_IC+p.E1_iNa_IC_noise*rand(1,p.E1_Npop);
E1_iNa_h = zeros(nsamp,p.E1_Npop);
E1_iNa_h(1,:) = E1_iNa_h_last;
E1_iK_n_last = p.E1_iK_n_IC+p.E1_iK_IC_noise*rand(1,p.E1_Npop);
E1_iK_n = zeros(nsamp,p.E1_Npop);
E1_iK_n(1,:) = E1_iK_n_last;
I1_V_last =  -rand(1, p.I1_Npop)*20;
I1_V = zeros(nsamp,p.I1_Npop);
I1_V(1,:) = I1_V_last;
I1_iNa_m_last = p.I1_iNa_m_IC+p.I1_iNa_IC_noise*rand(1,p.I1_Npop);
I1_iNa_m = zeros(nsamp,p.I1_Npop);
I1_iNa_m(1,:) = I1_iNa_m_last;
I1_iNa_h_last = p.I1_iNa_h_IC+p.I1_iNa_IC_noise*rand(1,p.I1_Npop);
I1_iNa_h = zeros(nsamp,p.I1_Npop);
I1_iNa_h(1,:) = I1_iNa_h_last;
I1_iK_n_last = p.I1_iK_n_IC+p.I1_iK_IC_noise*rand(1,p.I1_Npop);
I1_iK_n = zeros(nsamp,p.I1_Npop);
I1_iK_n(1,:) = I1_iK_n_last;
S1_V_last =  -rand(1, p.S1_Npop)*20;
S1_V = zeros(nsamp,p.S1_Npop);
S1_V(1,:) = S1_V_last;
S1_iNa_m_last = p.S1_iNa_m_IC+p.S1_iNa_IC_noise*rand(1,p.S1_Npop);
S1_iNa_m = zeros(nsamp,p.S1_Npop);
S1_iNa_m(1,:) = S1_iNa_m_last;
S1_iNa_h_last = p.S1_iNa_h_IC+p.S1_iNa_IC_noise*rand(1,p.S1_Npop);
S1_iNa_h = zeros(nsamp,p.S1_Npop);
S1_iNa_h(1,:) = S1_iNa_h_last;
S1_iK_n_last = p.S1_iK_n_IC+p.S1_iK_IC_noise*rand(1,p.S1_Npop);
S1_iK_n = zeros(nsamp,p.S1_Npop);
S1_iK_n(1,:) = S1_iK_n_last;
S2_V_last =  -rand(1, p.S2_Npop)*20;
S2_V = zeros(nsamp,p.S2_Npop);
S2_V(1,:) = S2_V_last;
S2_iNa_m_last = p.S2_iNa_m_IC+p.S2_iNa_IC_noise*rand(1,p.S2_Npop);
S2_iNa_m = zeros(nsamp,p.S2_Npop);
S2_iNa_m(1,:) = S2_iNa_m_last;
S2_iNa_h_last = p.S2_iNa_h_IC+p.S2_iNa_IC_noise*rand(1,p.S2_Npop);
S2_iNa_h = zeros(nsamp,p.S2_Npop);
S2_iNa_h(1,:) = S2_iNa_h_last;
S2_iK_n_last = p.S2_iK_n_IC+p.S2_iK_IC_noise*rand(1,p.S2_Npop);
S2_iK_n = zeros(nsamp,p.S2_Npop);
S2_iK_n(1,:) = S2_iK_n_last;
I1_E1_iAMPActx_s_last =  p.I1_E1_iAMPActx_IC+p.I1_E1_iAMPActx_IC_noise.*rand(1,p.E1_Npop);
I1_E1_iAMPActx_s = zeros(nsamp,p.E1_Npop);
I1_E1_iAMPActx_s(1,:) = I1_E1_iAMPActx_s_last;
E1_E1_iAMPActx_s_last =  p.E1_E1_iAMPActx_IC+p.E1_E1_iAMPActx_IC_noise.*rand(1,p.E1_Npop);
E1_E1_iAMPActx_s = zeros(nsamp,p.E1_Npop);
E1_E1_iAMPActx_s(1,:) = E1_E1_iAMPActx_s_last;
E1_I1_iGABActx_s_last =  p.E1_I1_iGABActx_IC+p.E1_I1_iGABActx_IC_noise.*rand(1,p.I1_Npop);
E1_I1_iGABActx_s = zeros(nsamp,p.I1_Npop);
E1_I1_iGABActx_s(1,:) = E1_I1_iGABActx_s_last;
I1_I1_iGABActx_s_last =  p.I1_I1_iGABActx_IC+p.I1_I1_iGABActx_IC_noise.*rand(1,p.I1_Npop);
I1_I1_iGABActx_s = zeros(nsamp,p.I1_Npop);
I1_I1_iGABActx_s(1,:) = I1_I1_iGABActx_s_last;
S2_S1_iAMPActx_s_last =  p.S2_S1_iAMPActx_IC+p.S2_S1_iAMPActx_IC_noise.*rand(1,p.S1_Npop);
S2_S1_iAMPActx_s = zeros(nsamp,p.S1_Npop);
S2_S1_iAMPActx_s(1,:) = S2_S1_iAMPActx_s_last;
S1_S1_iAMPActx_s_last =  p.S1_S1_iAMPActx_IC+p.S1_S1_iAMPActx_IC_noise.*rand(1,p.S1_Npop);
S1_S1_iAMPActx_s = zeros(nsamp,p.S1_Npop);
S1_S1_iAMPActx_s(1,:) = S1_S1_iAMPActx_s_last;
S1_S2_iAMPActx_s_last =  p.S1_S2_iAMPActx_IC+p.S1_S2_iAMPActx_IC_noise.*rand(1,p.S2_Npop);
S1_S2_iAMPActx_s = zeros(nsamp,p.S2_Npop);
S1_S2_iAMPActx_s(1,:) = S1_S2_iAMPActx_s_last;
S2_S2_iAMPActx_s_last =  p.S2_S2_iAMPActx_IC+p.S2_S2_iAMPActx_IC_noise.*rand(1,p.S2_Npop);
S2_S2_iAMPActx_s = zeros(nsamp,p.S2_Npop);
S2_S2_iAMPActx_s(1,:) = S2_S2_iAMPActx_s_last;
E1_S1_iAMPActx_s_last =  p.E1_S1_iAMPActx_IC+p.E1_S1_iAMPActx_IC_noise.*rand(1,p.S1_Npop);
E1_S1_iAMPActx_s = zeros(nsamp,p.S1_Npop);
E1_S1_iAMPActx_s(1,:) = E1_S1_iAMPActx_s_last;
E1_S2_iAMPActx_s_last =  p.E1_S2_iAMPActx_IC+p.E1_S2_iAMPActx_IC_noise.*rand(1,p.S2_Npop);
E1_S2_iAMPActx_s = zeros(nsamp,p.S2_Npop);
E1_S2_iAMPActx_s(1,:) = E1_S2_iAMPActx_s_last;

% MONITORS:
E1_ctx_iPoisson_g_poisson = zeros(nsamp,p.E1_Npop);
E1_ctx_iPoisson_g_poisson(1,:)=p.E1_ctx_iPoisson_g_poisson.*E1_ctx_iPoisson_s_poisson(1+floor(t/(0.5*dt)),:);
E1_ctx_iPoisson_I_poisson = zeros(nsamp,p.E1_Npop);
E1_ctx_iPoisson_I_poisson(1,:)=-((p.E1_ctx_iPoisson_g_poisson.*E1_ctx_iPoisson_s_poisson(1+floor(t/(0.5*dt)),:))).*(E1_V_last-p.E1_ctx_iPoisson_E_poisson);
I1_ctx_iPoisson_g_poisson = zeros(nsamp,p.I1_Npop);
I1_ctx_iPoisson_g_poisson(1,:)=p.I1_ctx_iPoisson_g_poisson.*I1_ctx_iPoisson_s_poisson(1+floor(t/(0.5*dt)),:);
I1_ctx_iPoisson_I_poisson = zeros(nsamp,p.I1_Npop);
I1_ctx_iPoisson_I_poisson(1,:)=-((p.I1_ctx_iPoisson_g_poisson.*I1_ctx_iPoisson_s_poisson(1+floor(t/(0.5*dt)),:))).*(I1_V_last-p.I1_ctx_iPoisson_E_poisson);
S1_ctx_iPoisson_g_poisson = zeros(nsamp,p.S1_Npop);
S1_ctx_iPoisson_g_poisson(1,:)=p.S1_ctx_iPoisson_g_poisson.*S1_ctx_iPoisson_s_poisson(1+floor(t/(0.5*dt)),:);
S1_ctx_iPoisson_I_poisson = zeros(nsamp,p.S1_Npop);
S1_ctx_iPoisson_I_poisson(1,:)=-((p.S1_ctx_iPoisson_g_poisson.*S1_ctx_iPoisson_s_poisson(1+floor(t/(0.5*dt)),:))).*(S1_V_last-p.S1_ctx_iPoisson_E_poisson);
S2_ctx_iPoisson_g_poisson = zeros(nsamp,p.S2_Npop);
S2_ctx_iPoisson_g_poisson(1,:)=p.S2_ctx_iPoisson_g_poisson.*S2_ctx_iPoisson_s_poisson(1+floor(t/(0.5*dt)),:);
S2_ctx_iPoisson_I_poisson = zeros(nsamp,p.S2_Npop);
S2_ctx_iPoisson_I_poisson(1,:)=-((p.S2_ctx_iPoisson_g_poisson.*S2_ctx_iPoisson_s_poisson(1+floor(t/(0.5*dt)),:))).*(S2_V_last-p.S2_ctx_iPoisson_E_poisson);

% ###########################################################
% Numerical integration:
% ###########################################################
n=2;
for k=2:ntime
  t=T(k-1);
  E1_V_k1 = (p.E1_Iapp + ((((-p.E1_iNa_gNa.*E1_iNa_m_last.^3.*E1_iNa_h_last.*(E1_V_last-p.E1_iNa_ENa))))+((((-p.E1_iK_gK.*E1_iK_n_last.^4.*(E1_V_last-p.E1_iK_EK))))+((((-((p.E1_ctx_iPoisson_g_poisson.*E1_ctx_iPoisson_s_poisson(1+floor(t/(0.5*dt)),:))).*(E1_V_last-p.E1_ctx_iPoisson_E_poisson))))+((-(((p.E1_E1_iAMPActx_gAMPA.*(E1_E1_iAMPActx_s_last*p.E1_E1_iAMPActx_netcon).*(E1_V_last-p.E1_E1_iAMPActx_EAMPA)))))+((((-p.E1_I1_iGABActx_gGABAa.*(E1_I1_iGABActx_s_last*p.E1_I1_iGABActx_netcon).*(E1_V_last-p.E1_I1_iGABActx_EGABAa))))+((-(((p.E1_S1_iAMPActx_gAMPA.*(E1_S1_iAMPActx_s_last*p.E1_S1_iAMPActx_netcon).*(E1_V_last-p.E1_S1_iAMPActx_EAMPA)))))+((-(((p.E1_S2_iAMPActx_gAMPA.*(E1_S2_iAMPActx_s_last*p.E1_S2_iAMPActx_netcon).*(E1_V_last-p.E1_S2_iAMPActx_EAMPA)))))))))))) + p.E1_noise*randn(1, p.E1_Npop))/p.E1_C;
  E1_iNa_m_k1 = (((2.5-.1*(E1_V_last+65))./(exp(2.5-.1*(E1_V_last+65))-1))).*(1-E1_iNa_m_last)-((4*exp(-(E1_V_last+65)/18))).*E1_iNa_m_last;
  E1_iNa_h_k1 = ((.07*exp(-(E1_V_last+65)/20))).*(1-E1_iNa_h_last)-((1./(exp(3-.1*(E1_V_last+65))+1))).*E1_iNa_h_last;
  E1_iK_n_k1 = (((.1-.01*(E1_V_last+65))./(exp(1-.1*(E1_V_last+65))-1))).*(1-E1_iK_n_last)-((.125*exp(-(E1_V_last+65)/80))).*E1_iK_n_last;
  I1_V_k1 = (p.I1_Iapp + ((((-p.I1_iNa_gNa.*I1_iNa_m_last.^3.*I1_iNa_h_last.*(I1_V_last-p.I1_iNa_ENa))))+((((-p.I1_iK_gK.*I1_iK_n_last.^4.*(I1_V_last-p.I1_iK_EK))))+((((-((p.I1_ctx_iPoisson_g_poisson.*I1_ctx_iPoisson_s_poisson(1+floor(t/(0.5*dt)),:))).*(I1_V_last-p.I1_ctx_iPoisson_E_poisson))))+((-(((p.I1_E1_iAMPActx_gAMPA.*(I1_E1_iAMPActx_s_last*p.I1_E1_iAMPActx_netcon).*(I1_V_last-p.I1_E1_iAMPActx_EAMPA)))))+((((-p.I1_I1_iGABActx_gGABAa.*(I1_I1_iGABActx_s_last*p.I1_I1_iGABActx_netcon).*(I1_V_last-p.I1_I1_iGABActx_EGABAa))))))))) + p.I1_noise*randn(1, p.I1_Npop))/p.I1_C;
  I1_iNa_m_k1 = (((2.5-.1*(I1_V_last+65))./(exp(2.5-.1*(I1_V_last+65))-1))).*(1-I1_iNa_m_last)-((4*exp(-(I1_V_last+65)/18))).*I1_iNa_m_last;
  I1_iNa_h_k1 = ((.07*exp(-(I1_V_last+65)/20))).*(1-I1_iNa_h_last)-((1./(exp(3-.1*(I1_V_last+65))+1))).*I1_iNa_h_last;
  I1_iK_n_k1 = (((.1-.01*(I1_V_last+65))./(exp(1-.1*(I1_V_last+65))-1))).*(1-I1_iK_n_last)-((.125*exp(-(I1_V_last+65)/80))).*I1_iK_n_last;
  S1_V_k1 = (p.S1_Iapp + ((((-p.S1_iNa_gNa.*S1_iNa_m_last.^3.*S1_iNa_h_last.*(S1_V_last-p.S1_iNa_ENa))))+((((-p.S1_iK_gK.*S1_iK_n_last.^4.*(S1_V_last-p.S1_iK_EK))))+((((-((p.S1_ctx_iPoisson_g_poisson.*S1_ctx_iPoisson_s_poisson(1+floor(t/(0.5*dt)),:))).*(S1_V_last-p.S1_ctx_iPoisson_E_poisson))))+((-(((p.S1_S1_iAMPActx_gAMPA.*(S1_S1_iAMPActx_s_last*p.S1_S1_iAMPActx_netcon).*(S1_V_last-p.S1_S1_iAMPActx_EAMPA)))))+((-(((p.S1_S2_iAMPActx_gAMPA.*(S1_S2_iAMPActx_s_last*p.S1_S2_iAMPActx_netcon).*(S1_V_last-p.S1_S2_iAMPActx_EAMPA)))))))))) + p.S1_noise*randn(1, p.S1_Npop))/p.S1_C;
  S1_iNa_m_k1 = (((2.5-.1*(S1_V_last+65))./(exp(2.5-.1*(S1_V_last+65))-1))).*(1-S1_iNa_m_last)-((4*exp(-(S1_V_last+65)/18))).*S1_iNa_m_last;
  S1_iNa_h_k1 = ((.07*exp(-(S1_V_last+65)/20))).*(1-S1_iNa_h_last)-((1./(exp(3-.1*(S1_V_last+65))+1))).*S1_iNa_h_last;
  S1_iK_n_k1 = (((.1-.01*(S1_V_last+65))./(exp(1-.1*(S1_V_last+65))-1))).*(1-S1_iK_n_last)-((.125*exp(-(S1_V_last+65)/80))).*S1_iK_n_last;
  S2_V_k1 = (p.S2_Iapp + ((((-p.S2_iNa_gNa.*S2_iNa_m_last.^3.*S2_iNa_h_last.*(S2_V_last-p.S2_iNa_ENa))))+((((-p.S2_iK_gK.*S2_iK_n_last.^4.*(S2_V_last-p.S2_iK_EK))))+((((-((p.S2_ctx_iPoisson_g_poisson.*S2_ctx_iPoisson_s_poisson(1+floor(t/(0.5*dt)),:))).*(S2_V_last-p.S2_ctx_iPoisson_E_poisson))))+((-(((p.S2_S1_iAMPActx_gAMPA.*(S2_S1_iAMPActx_s_last*p.S2_S1_iAMPActx_netcon).*(S2_V_last-p.S2_S1_iAMPActx_EAMPA)))))+((-(((p.S2_S2_iAMPActx_gAMPA.*(S2_S2_iAMPActx_s_last*p.S2_S2_iAMPActx_netcon).*(S2_V_last-p.S2_S2_iAMPActx_EAMPA)))))))))) + p.S2_noise*randn(1, p.S2_Npop))/p.S2_C;
  S2_iNa_m_k1 = (((2.5-.1*(S2_V_last+65))./(exp(2.5-.1*(S2_V_last+65))-1))).*(1-S2_iNa_m_last)-((4*exp(-(S2_V_last+65)/18))).*S2_iNa_m_last;
  S2_iNa_h_k1 = ((.07*exp(-(S2_V_last+65)/20))).*(1-S2_iNa_h_last)-((1./(exp(3-.1*(S2_V_last+65))+1))).*S2_iNa_h_last;
  S2_iK_n_k1 = (((.1-.01*(S2_V_last+65))./(exp(1-.1*(S2_V_last+65))-1))).*(1-S2_iK_n_last)-((.125*exp(-(S2_V_last+65)/80))).*S2_iK_n_last;
  I1_E1_iAMPActx_s_k1 = -I1_E1_iAMPActx_s_last./p.I1_E1_iAMPActx_tauAMPA + ((1-I1_E1_iAMPActx_s_last)/p.I1_E1_iAMPActx_tauAMPAr).*(1+tanh(E1_V_last/10));
  E1_E1_iAMPActx_s_k1 = -E1_E1_iAMPActx_s_last./p.E1_E1_iAMPActx_tauAMPA + ((1-E1_E1_iAMPActx_s_last)/p.E1_E1_iAMPActx_tauAMPAr).*(1+tanh(E1_V_last/10));
  E1_I1_iGABActx_s_k1 = -E1_I1_iGABActx_s_last./p.E1_I1_iGABActx_tauGABA + 1/2*(1+tanh(I1_V_last/10)).*((1-E1_I1_iGABActx_s_last)/p.E1_I1_iGABActx_tauGABAr);
  I1_I1_iGABActx_s_k1 = -I1_I1_iGABActx_s_last./p.I1_I1_iGABActx_tauGABA + 1/2*(1+tanh(I1_V_last/10)).*((1-I1_I1_iGABActx_s_last)/p.I1_I1_iGABActx_tauGABAr);
  S2_S1_iAMPActx_s_k1 = -S2_S1_iAMPActx_s_last./p.S2_S1_iAMPActx_tauAMPA + ((1-S2_S1_iAMPActx_s_last)/p.S2_S1_iAMPActx_tauAMPAr).*(1+tanh(S1_V_last/10));
  S1_S1_iAMPActx_s_k1 = -S1_S1_iAMPActx_s_last./p.S1_S1_iAMPActx_tauAMPA + ((1-S1_S1_iAMPActx_s_last)/p.S1_S1_iAMPActx_tauAMPAr).*(1+tanh(S1_V_last/10));
  S1_S2_iAMPActx_s_k1 = -S1_S2_iAMPActx_s_last./p.S1_S2_iAMPActx_tauAMPA + ((1-S1_S2_iAMPActx_s_last)/p.S1_S2_iAMPActx_tauAMPAr).*(1+tanh(S2_V_last/10));
  S2_S2_iAMPActx_s_k1 = -S2_S2_iAMPActx_s_last./p.S2_S2_iAMPActx_tauAMPA + ((1-S2_S2_iAMPActx_s_last)/p.S2_S2_iAMPActx_tauAMPAr).*(1+tanh(S2_V_last/10));
  E1_S1_iAMPActx_s_k1 = -E1_S1_iAMPActx_s_last./p.E1_S1_iAMPActx_tauAMPA + ((1-E1_S1_iAMPActx_s_last)/p.E1_S1_iAMPActx_tauAMPAr).*(1+tanh(S1_V_last/10));
  E1_S2_iAMPActx_s_k1 = -E1_S2_iAMPActx_s_last./p.E1_S2_iAMPActx_tauAMPA + ((1-E1_S2_iAMPActx_s_last)/p.E1_S2_iAMPActx_tauAMPAr).*(1+tanh(S2_V_last/10));

  % ------------------------------------------------------------
  % Update state variables:
  % ------------------------------------------------------------
  E1_V_last = E1_V_last+dt*E1_V_k1;
  E1_iNa_m_last = E1_iNa_m_last+dt*E1_iNa_m_k1;
  E1_iNa_h_last = E1_iNa_h_last+dt*E1_iNa_h_k1;
  E1_iK_n_last = E1_iK_n_last+dt*E1_iK_n_k1;
  I1_V_last = I1_V_last+dt*I1_V_k1;
  I1_iNa_m_last = I1_iNa_m_last+dt*I1_iNa_m_k1;
  I1_iNa_h_last = I1_iNa_h_last+dt*I1_iNa_h_k1;
  I1_iK_n_last = I1_iK_n_last+dt*I1_iK_n_k1;
  S1_V_last = S1_V_last+dt*S1_V_k1;
  S1_iNa_m_last = S1_iNa_m_last+dt*S1_iNa_m_k1;
  S1_iNa_h_last = S1_iNa_h_last+dt*S1_iNa_h_k1;
  S1_iK_n_last = S1_iK_n_last+dt*S1_iK_n_k1;
  S2_V_last = S2_V_last+dt*S2_V_k1;
  S2_iNa_m_last = S2_iNa_m_last+dt*S2_iNa_m_k1;
  S2_iNa_h_last = S2_iNa_h_last+dt*S2_iNa_h_k1;
  S2_iK_n_last = S2_iK_n_last+dt*S2_iK_n_k1;
  I1_E1_iAMPActx_s_last = I1_E1_iAMPActx_s_last+dt*I1_E1_iAMPActx_s_k1;
  E1_E1_iAMPActx_s_last = E1_E1_iAMPActx_s_last+dt*E1_E1_iAMPActx_s_k1;
  E1_I1_iGABActx_s_last = E1_I1_iGABActx_s_last+dt*E1_I1_iGABActx_s_k1;
  I1_I1_iGABActx_s_last = I1_I1_iGABActx_s_last+dt*I1_I1_iGABActx_s_k1;
  S2_S1_iAMPActx_s_last = S2_S1_iAMPActx_s_last+dt*S2_S1_iAMPActx_s_k1;
  S1_S1_iAMPActx_s_last = S1_S1_iAMPActx_s_last+dt*S1_S1_iAMPActx_s_k1;
  S1_S2_iAMPActx_s_last = S1_S2_iAMPActx_s_last+dt*S1_S2_iAMPActx_s_k1;
  S2_S2_iAMPActx_s_last = S2_S2_iAMPActx_s_last+dt*S2_S2_iAMPActx_s_k1;
  E1_S1_iAMPActx_s_last = E1_S1_iAMPActx_s_last+dt*E1_S1_iAMPActx_s_k1;
  E1_S2_iAMPActx_s_last = E1_S2_iAMPActx_s_last+dt*E1_S2_iAMPActx_s_k1;

  if mod(k,downsample_factor)==0 % store this time point

  % ------------------------------------------------------------
  % Store state variables:
  % ------------------------------------------------------------
    E1_V(n,:) = E1_V_last;
    E1_iNa_m(n,:) = E1_iNa_m_last;
    E1_iNa_h(n,:) = E1_iNa_h_last;
    E1_iK_n(n,:) = E1_iK_n_last;
    I1_V(n,:) = I1_V_last;
    I1_iNa_m(n,:) = I1_iNa_m_last;
    I1_iNa_h(n,:) = I1_iNa_h_last;
    I1_iK_n(n,:) = I1_iK_n_last;
    S1_V(n,:) = S1_V_last;
    S1_iNa_m(n,:) = S1_iNa_m_last;
    S1_iNa_h(n,:) = S1_iNa_h_last;
    S1_iK_n(n,:) = S1_iK_n_last;
    S2_V(n,:) = S2_V_last;
    S2_iNa_m(n,:) = S2_iNa_m_last;
    S2_iNa_h(n,:) = S2_iNa_h_last;
    S2_iK_n(n,:) = S2_iK_n_last;
    I1_E1_iAMPActx_s(n,:) = I1_E1_iAMPActx_s_last;
    E1_E1_iAMPActx_s(n,:) = E1_E1_iAMPActx_s_last;
    E1_I1_iGABActx_s(n,:) = E1_I1_iGABActx_s_last;
    I1_I1_iGABActx_s(n,:) = I1_I1_iGABActx_s_last;
    S2_S1_iAMPActx_s(n,:) = S2_S1_iAMPActx_s_last;
    S1_S1_iAMPActx_s(n,:) = S1_S1_iAMPActx_s_last;
    S1_S2_iAMPActx_s(n,:) = S1_S2_iAMPActx_s_last;
    S2_S2_iAMPActx_s(n,:) = S2_S2_iAMPActx_s_last;
    E1_S1_iAMPActx_s(n,:) = E1_S1_iAMPActx_s_last;
    E1_S2_iAMPActx_s(n,:) = E1_S2_iAMPActx_s_last;

  % ------------------------------------------------------------
  % Update monitors:
  % ------------------------------------------------------------
  E1_ctx_iPoisson_g_poisson(n,:)=p.E1_ctx_iPoisson_g_poisson.*E1_ctx_iPoisson_s_poisson(1+floor(t/(0.5*dt)),:);
  E1_ctx_iPoisson_I_poisson(n,:)=-((p.E1_ctx_iPoisson_g_poisson.*E1_ctx_iPoisson_s_poisson(1+floor(t/(0.5*dt)),:))).*(E1_V(n,:)-p.E1_ctx_iPoisson_E_poisson);
  I1_ctx_iPoisson_g_poisson(n,:)=p.I1_ctx_iPoisson_g_poisson.*I1_ctx_iPoisson_s_poisson(1+floor(t/(0.5*dt)),:);
  I1_ctx_iPoisson_I_poisson(n,:)=-((p.I1_ctx_iPoisson_g_poisson.*I1_ctx_iPoisson_s_poisson(1+floor(t/(0.5*dt)),:))).*(I1_V(n,:)-p.I1_ctx_iPoisson_E_poisson);
  S1_ctx_iPoisson_g_poisson(n,:)=p.S1_ctx_iPoisson_g_poisson.*S1_ctx_iPoisson_s_poisson(1+floor(t/(0.5*dt)),:);
  S1_ctx_iPoisson_I_poisson(n,:)=-((p.S1_ctx_iPoisson_g_poisson.*S1_ctx_iPoisson_s_poisson(1+floor(t/(0.5*dt)),:))).*(S1_V(n,:)-p.S1_ctx_iPoisson_E_poisson);
  S2_ctx_iPoisson_g_poisson(n,:)=p.S2_ctx_iPoisson_g_poisson.*S2_ctx_iPoisson_s_poisson(1+floor(t/(0.5*dt)),:);
  S2_ctx_iPoisson_I_poisson(n,:)=-((p.S2_ctx_iPoisson_g_poisson.*S2_ctx_iPoisson_s_poisson(1+floor(t/(0.5*dt)),:))).*(S2_V(n,:)-p.S2_ctx_iPoisson_E_poisson);

    n=n+1;
  end
end

T=T(1:downsample_factor:ntime);

end
