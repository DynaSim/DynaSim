function [T,E_V,E_iNa_m_na,E_iNa_h_na,E_iK_m_k,I_V,I_iNa_m_na,I_iNa_h_na,I_iK_m_k,E_E_iAMPA_s,I_E_iAMPA_s,E_I_iGABAa_s,I_I_iGABAa_s,E_current,E_iE,E_iI,E_iL_i_l,E_iNa_i_Na,E_iK_i_k,E_iACh_i_ACh,E_iNE_i_NE,I_current,I_iE,I_iI,I_iL_i_l,I_iNa_i_Na,I_iK_i_k,E_iPoisson_s,E_iACh_eta_ACh,E_iACh_rnd_sample,E_iACh_M_ACh,E_iNE_eta_NE,E_iNE_rnd_sample,E_iNE_M_NE,I_iPoisson_s]=solve_ode

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

% the 'shuffle' random seed has been set in advance.

% ------------------------------------------------------------
% Fixed variables:
% ------------------------------------------------------------
E_iPoisson_s = getPoissonGating(p.E_iPoisson_baseline,p.E_iPoisson_DC,p.E_iPoisson_AC,p.E_iPoisson_f,p.E_iPoisson_phi,p.E_iPoisson_onset,p.E_iPoisson_offset,p.E_iPoisson_tau,T,p.E_Npop);
E_iACh_eta_ACh =  rand(1, p.E_Npop) <= p.E_iACh_rho_ACh;
E_iACh_rnd_sample =  randn(1, p.E_Npop);
E_iACh_M_ACh =  p.E_iACh_A_ACh*(1 + E_iACh_rnd_sample).*(E_iACh_rnd_sample > -1);
E_iNE_eta_NE =  rand(1, p.E_Npop) <= p.E_iNE_rho_NE;
E_iNE_rnd_sample =  randn(1, p.E_Npop);
E_iNE_M_NE =  p.E_iNE_A_NE*(1 + E_iNE_rnd_sample).*(E_iNE_rnd_sample > -1);
I_iPoisson_s = getPoissonGating(p.I_iPoisson_baseline,p.I_iPoisson_DC,p.I_iPoisson_AC,p.I_iPoisson_f,p.I_iPoisson_phi,p.I_iPoisson_onset,p.I_iPoisson_offset,p.I_iPoisson_tau,T,p.I_Npop);

% ------------------------------------------------------------
% Initial conditions:
% ------------------------------------------------------------
t=0; k=1;

% STATE_VARIABLES:
E_V_last =  -90 + 20*randn(1, p.E_Npop);
E_V = zeros(nsamp,p.E_Npop);
E_V(1,:) = E_V_last;
E_iNa_m_na_last =  0 * ones(1,p.E_Npop);
E_iNa_m_na = zeros(nsamp,p.E_Npop);
E_iNa_m_na(1,:) = E_iNa_m_na_last;
E_iNa_h_na_last =  0 * ones(1,p.E_Npop);
E_iNa_h_na = zeros(nsamp,p.E_Npop);
E_iNa_h_na(1,:) = E_iNa_h_na_last;
E_iK_m_k_last =  0 * ones(1,p.E_Npop);
E_iK_m_k = zeros(nsamp,p.E_Npop);
E_iK_m_k(1,:) = E_iK_m_k_last;
I_V_last =  -90 + 20*randn(1, p.I_Npop);
I_V = zeros(nsamp,p.I_Npop);
I_V(1,:) = I_V_last;
I_iNa_m_na_last =  0 * ones(1,p.I_Npop);
I_iNa_m_na = zeros(nsamp,p.I_Npop);
I_iNa_m_na(1,:) = I_iNa_m_na_last;
I_iNa_h_na_last =  0 * ones(1,p.I_Npop);
I_iNa_h_na = zeros(nsamp,p.I_Npop);
I_iNa_h_na(1,:) = I_iNa_h_na_last;
I_iK_m_k_last =  0 * ones(1,p.I_Npop);
I_iK_m_k = zeros(nsamp,p.I_Npop);
I_iK_m_k(1,:) = I_iK_m_k_last;
E_E_iAMPA_s_last =  p.E_E_iAMPA_IC+p.E_E_iAMPA_IC_noise.*rand(1,p.E_Npop);
E_E_iAMPA_s = zeros(nsamp,p.E_Npop);
E_E_iAMPA_s(1,:) = E_E_iAMPA_s_last;
I_E_iAMPA_s_last =  p.I_E_iAMPA_IC+p.I_E_iAMPA_IC_noise.*rand(1,p.E_Npop);
I_E_iAMPA_s = zeros(nsamp,p.E_Npop);
I_E_iAMPA_s(1,:) = I_E_iAMPA_s_last;
E_I_iGABAa_s_last =  p.E_I_iGABAa_IC+p.E_I_iGABAa_IC_noise.*rand(1,p.I_Npop);
E_I_iGABAa_s = zeros(nsamp,p.I_Npop);
E_I_iGABAa_s(1,:) = E_I_iGABAa_s_last;
I_I_iGABAa_s_last =  p.I_I_iGABAa_IC+p.I_I_iGABAa_IC_noise.*rand(1,p.I_Npop);
I_I_iGABAa_s = zeros(nsamp,p.I_Npop);
I_I_iGABAa_s(1,:) = I_I_iGABAa_s_last;

% MONITORS:
E_current = zeros(nsamp,p.E_Npop);
E_current(1,:)=((((-p.E_iL_g_l*(E_V_last-p.E_iL_E_l))))+((((-p.E_iNa_g_na*E_iNa_m_na_last.^3.*E_iNa_h_na_last.*(E_V_last-p.E_iNa_E_na))))+((((-p.E_iK_g_k*E_iK_m_k_last.^4.*(E_V_last-p.E_iK_E_k))))+((((-((p.E_iPoisson_gext.*E_iPoisson_s(k,:))).*(E_V_last-p.E_iPoisson_Eext))))+((((p.E_iACh_Y_ACh*E_iACh_eta_ACh.*E_iACh_M_ACh)))+((((E_iNE_eta_NE.*E_iNE_M_NE)))+((((-p.E_E_iAMPA_gAMPA.*(E_E_iAMPA_s_last*p.E_E_iAMPA_netcon).*(E_V_last-p.E_E_iAMPA_EAMPA))))+((((-p.E_I_iGABAa_gGABAa.*(E_I_iGABAa_s_last*p.E_I_iGABAa_netcon).*(E_V_last-p.E_I_iGABAa_EGABAa))))))))))));
E_iE = zeros(nsamp,p.E_Npop);
E_iE(1,:)=((((-((p.E_iPoisson_gext.*E_iPoisson_s(k,:))).*(E_V_last-p.E_iPoisson_Eext))))+((((-p.E_E_iAMPA_gAMPA.*(E_E_iAMPA_s_last*p.E_E_iAMPA_netcon).*(E_V_last-p.E_E_iAMPA_EAMPA))))));
E_iI = zeros(nsamp,p.E_Npop);
E_iI(1,:)=((((-p.E_I_iGABAa_gGABAa.*(E_I_iGABAa_s_last*p.E_I_iGABAa_netcon).*(E_V_last-p.E_I_iGABAa_EGABAa)))));
E_iL_i_l = zeros(nsamp,p.E_Npop);
E_iL_i_l(1,:)=-p.E_iL_g_l*(E_V_last-p.E_iL_E_l);
E_iNa_i_Na = zeros(nsamp,p.E_Npop);
E_iNa_i_Na(1,:)=-p.E_iNa_g_na*E_iNa_m_na_last.^3.*E_iNa_h_na_last.*(E_V_last-p.E_iNa_E_na);
E_iK_i_k = zeros(nsamp,p.E_Npop);
E_iK_i_k(1,:)=-p.E_iK_g_k*E_iK_m_k_last.^4.*(E_V_last-p.E_iK_E_k);
E_iACh_i_ACh = zeros(nsamp,p.E_Npop);
E_iACh_i_ACh(1,:)=p.E_iACh_Y_ACh*E_iACh_eta_ACh.*E_iACh_M_ACh;
E_iNE_i_NE = zeros(nsamp,p.E_Npop);
E_iNE_i_NE(1,:)=E_iNE_eta_NE.*E_iNE_M_NE;
I_current = zeros(nsamp,p.I_Npop);
I_current(1,:)=((((-p.I_iL_g_l*(I_V_last-p.I_iL_E_l))))+((((-p.I_iNa_g_na*I_iNa_m_na_last.^3.*I_iNa_h_na_last.*(I_V_last-p.I_iNa_E_na))))+((((-p.I_iK_g_k*I_iK_m_k_last.^4.*(I_V_last-p.I_iK_E_k))))+((((-((p.I_iPoisson_gext.*I_iPoisson_s(k,:))).*(I_V_last-p.I_iPoisson_Eext))))+((((-p.I_E_iAMPA_gAMPA.*(I_E_iAMPA_s_last*p.I_E_iAMPA_netcon).*(I_V_last-p.I_E_iAMPA_EAMPA))))+((((-p.I_I_iGABAa_gGABAa.*(I_I_iGABAa_s_last*p.I_I_iGABAa_netcon).*(I_V_last-p.I_I_iGABAa_EGABAa))))))))));
I_iE = zeros(nsamp,p.I_Npop);
I_iE(1,:)=((((-((p.I_iPoisson_gext.*I_iPoisson_s(k,:))).*(I_V_last-p.I_iPoisson_Eext))))+((((-p.I_E_iAMPA_gAMPA.*(I_E_iAMPA_s_last*p.I_E_iAMPA_netcon).*(I_V_last-p.I_E_iAMPA_EAMPA))))));
I_iI = zeros(nsamp,p.I_Npop);
I_iI(1,:)=((((-p.I_I_iGABAa_gGABAa.*(I_I_iGABAa_s_last*p.I_I_iGABAa_netcon).*(I_V_last-p.I_I_iGABAa_EGABAa)))));
I_iL_i_l = zeros(nsamp,p.I_Npop);
I_iL_i_l(1,:)=-p.I_iL_g_l*(I_V_last-p.I_iL_E_l);
I_iNa_i_Na = zeros(nsamp,p.I_Npop);
I_iNa_i_Na(1,:)=-p.I_iNa_g_na*I_iNa_m_na_last.^3.*I_iNa_h_na_last.*(I_V_last-p.I_iNa_E_na);
I_iK_i_k = zeros(nsamp,p.I_Npop);
I_iK_i_k(1,:)=-p.I_iK_g_k*I_iK_m_k_last.^4.*(I_V_last-p.I_iK_E_k);

% ###########################################################
% Memory check:
% ###########################################################
try 
  memoryUsed = memoryUsageCallerGB(); 
  fprintf('Total Memory Used <= %i GB \n', ceil(memoryUsed)); 
end 

% ###########################################################
% Numerical integration:
% ###########################################################
n=2;
for k=2:ntime
  t=T(k-1);
  E_V_k1 = ((((-p.E_iL_g_l*(E_V_last-p.E_iL_E_l))))+((((-p.E_iNa_g_na*E_iNa_m_na_last.^3.*E_iNa_h_na_last.*(E_V_last-p.E_iNa_E_na))))+((((-p.E_iK_g_k*E_iK_m_k_last.^4.*(E_V_last-p.E_iK_E_k))))+((((-((p.E_iPoisson_gext.*E_iPoisson_s(k,:))).*(E_V_last-p.E_iPoisson_Eext))))+((((p.E_iACh_Y_ACh*E_iACh_eta_ACh.*E_iACh_M_ACh)))+((((E_iNE_eta_NE.*E_iNE_M_NE)))+((((-p.E_E_iAMPA_gAMPA.*(E_E_iAMPA_s_last*p.E_E_iAMPA_netcon).*(E_V_last-p.E_E_iAMPA_EAMPA))))+((((-p.E_I_iGABAa_gGABAa.*(E_I_iGABAa_s_last*p.E_I_iGABAa_netcon).*(E_V_last-p.E_I_iGABAa_EGABAa))))))))))))/p.E_Cm;
  E_iNa_m_na_k1 = ((0.32*(E_V_last+54)./(1-exp(-(E_V_last+54)/4)))).*(1-E_iNa_m_na_last)-((0.28*(E_V_last+27)./(exp((E_V_last+27)/5)-1))).*E_iNa_m_na_last;
  E_iNa_h_na_k1 = ((0.128*exp(-(E_V_last+50)/18))).*(1-E_iNa_h_na_last)-((4./(1+exp(-(E_V_last+27)/5)))).*E_iNa_h_na_last;
  E_iK_m_k_k1 = ((0.032*(E_V_last+52)./(1-exp(-(E_V_last+52)/5)))).*(1-E_iK_m_k_last)-((0.5*exp(-(E_V_last+57)/40))).*E_iK_m_k_last;
  I_V_k1 = ((((-p.I_iL_g_l*(I_V_last-p.I_iL_E_l))))+((((-p.I_iNa_g_na*I_iNa_m_na_last.^3.*I_iNa_h_na_last.*(I_V_last-p.I_iNa_E_na))))+((((-p.I_iK_g_k*I_iK_m_k_last.^4.*(I_V_last-p.I_iK_E_k))))+((((-((p.I_iPoisson_gext.*I_iPoisson_s(k,:))).*(I_V_last-p.I_iPoisson_Eext))))+((((-p.I_E_iAMPA_gAMPA.*(I_E_iAMPA_s_last*p.I_E_iAMPA_netcon).*(I_V_last-p.I_E_iAMPA_EAMPA))))+((((-p.I_I_iGABAa_gGABAa.*(I_I_iGABAa_s_last*p.I_I_iGABAa_netcon).*(I_V_last-p.I_I_iGABAa_EGABAa))))))))))/p.I_Cm;
  I_iNa_m_na_k1 = ((0.32*(I_V_last+54)./(1-exp(-(I_V_last+54)/4)))).*(1-I_iNa_m_na_last)-((0.28*(I_V_last+27)./(exp((I_V_last+27)/5)-1))).*I_iNa_m_na_last;
  I_iNa_h_na_k1 = ((0.128*exp(-(I_V_last+50)/18))).*(1-I_iNa_h_na_last)-((4./(1+exp(-(I_V_last+27)/5)))).*I_iNa_h_na_last;
  I_iK_m_k_k1 = ((0.032*(I_V_last+52)./(1-exp(-(I_V_last+52)/5)))).*(1-I_iK_m_k_last)-((0.5*exp(-(I_V_last+57)/40))).*I_iK_m_k_last;
  E_E_iAMPA_s_k1 = -E_E_iAMPA_s_last./p.E_E_iAMPA_tauD + 1/2*(1+tanh(E_V_last/10)).*((1-E_E_iAMPA_s_last)/p.E_E_iAMPA_tauR);
  I_E_iAMPA_s_k1 = -I_E_iAMPA_s_last./p.I_E_iAMPA_tauD + 1/2*(1+tanh(E_V_last/10)).*((1-I_E_iAMPA_s_last)/p.I_E_iAMPA_tauR);
  E_I_iGABAa_s_k1 = -E_I_iGABAa_s_last./p.E_I_iGABAa_tauD + 1/2*(1+tanh(I_V_last/10)).*((1-E_I_iGABAa_s_last)/p.E_I_iGABAa_tauR);
  I_I_iGABAa_s_k1 = -I_I_iGABAa_s_last./p.I_I_iGABAa_tauD + 1/2*(1+tanh(I_V_last/10)).*((1-I_I_iGABAa_s_last)/p.I_I_iGABAa_tauR);

  % ------------------------------------------------------------
  % Update state variables:
  % ------------------------------------------------------------
  E_V_last = E_V_last+dt*E_V_k1;
  E_iNa_m_na_last = E_iNa_m_na_last+dt*E_iNa_m_na_k1;
  E_iNa_h_na_last = E_iNa_h_na_last+dt*E_iNa_h_na_k1;
  E_iK_m_k_last = E_iK_m_k_last+dt*E_iK_m_k_k1;
  I_V_last = I_V_last+dt*I_V_k1;
  I_iNa_m_na_last = I_iNa_m_na_last+dt*I_iNa_m_na_k1;
  I_iNa_h_na_last = I_iNa_h_na_last+dt*I_iNa_h_na_k1;
  I_iK_m_k_last = I_iK_m_k_last+dt*I_iK_m_k_k1;
  E_E_iAMPA_s_last = E_E_iAMPA_s_last+dt*E_E_iAMPA_s_k1;
  I_E_iAMPA_s_last = I_E_iAMPA_s_last+dt*I_E_iAMPA_s_k1;
  E_I_iGABAa_s_last = E_I_iGABAa_s_last+dt*E_I_iGABAa_s_k1;
  I_I_iGABAa_s_last = I_I_iGABAa_s_last+dt*I_I_iGABAa_s_k1;

  if mod(k,downsample_factor)==0 % store this time point

  % ------------------------------------------------------------
  % Store state variables:
  % ------------------------------------------------------------
    E_V(n,:) = E_V_last;
    E_iNa_m_na(n,:) = E_iNa_m_na_last;
    E_iNa_h_na(n,:) = E_iNa_h_na_last;
    E_iK_m_k(n,:) = E_iK_m_k_last;
    I_V(n,:) = I_V_last;
    I_iNa_m_na(n,:) = I_iNa_m_na_last;
    I_iNa_h_na(n,:) = I_iNa_h_na_last;
    I_iK_m_k(n,:) = I_iK_m_k_last;
    E_E_iAMPA_s(n,:) = E_E_iAMPA_s_last;
    I_E_iAMPA_s(n,:) = I_E_iAMPA_s_last;
    E_I_iGABAa_s(n,:) = E_I_iGABAa_s_last;
    I_I_iGABAa_s(n,:) = I_I_iGABAa_s_last;

  % ------------------------------------------------------------
  % Update monitors:
  % ------------------------------------------------------------
  E_current(n,:)=((((-p.E_iL_g_l*(E_V(n,:)-p.E_iL_E_l))))+((((-p.E_iNa_g_na*E_iNa_m_na(n,:).^3.*E_iNa_h_na(n,:).*(E_V(n,:)-p.E_iNa_E_na))))+((((-p.E_iK_g_k*E_iK_m_k(n,:).^4.*(E_V(n,:)-p.E_iK_E_k))))+((((-((p.E_iPoisson_gext.*E_iPoisson_s(k,:))).*(E_V(n,:)-p.E_iPoisson_Eext))))+((((p.E_iACh_Y_ACh*E_iACh_eta_ACh.*E_iACh_M_ACh)))+((((E_iNE_eta_NE.*E_iNE_M_NE)))+((((-p.E_E_iAMPA_gAMPA.*(E_E_iAMPA_s(n,:)*p.E_E_iAMPA_netcon).*(E_V(n,:)-p.E_E_iAMPA_EAMPA))))+((((-p.E_I_iGABAa_gGABAa.*(E_I_iGABAa_s(n,:)*p.E_I_iGABAa_netcon).*(E_V(n,:)-p.E_I_iGABAa_EGABAa))))))))))));
  E_iE(n,:)=((((-((p.E_iPoisson_gext.*E_iPoisson_s(k,:))).*(E_V(n,:)-p.E_iPoisson_Eext))))+((((-p.E_E_iAMPA_gAMPA.*(E_E_iAMPA_s(n,:)*p.E_E_iAMPA_netcon).*(E_V(n,:)-p.E_E_iAMPA_EAMPA))))));
  E_iI(n,:)=((((-p.E_I_iGABAa_gGABAa.*(E_I_iGABAa_s(n,:)*p.E_I_iGABAa_netcon).*(E_V(n,:)-p.E_I_iGABAa_EGABAa)))));
  E_iL_i_l(n,:)=-p.E_iL_g_l*(E_V(n,:)-p.E_iL_E_l);
  E_iNa_i_Na(n,:)=-p.E_iNa_g_na*E_iNa_m_na(n,:).^3.*E_iNa_h_na(n,:).*(E_V(n,:)-p.E_iNa_E_na);
  E_iK_i_k(n,:)=-p.E_iK_g_k*E_iK_m_k(n,:).^4.*(E_V(n,:)-p.E_iK_E_k);
  E_iACh_i_ACh(n,:)=p.E_iACh_Y_ACh*E_iACh_eta_ACh.*E_iACh_M_ACh;
  E_iNE_i_NE(n,:)=E_iNE_eta_NE.*E_iNE_M_NE;
  I_current(n,:)=((((-p.I_iL_g_l*(I_V(n,:)-p.I_iL_E_l))))+((((-p.I_iNa_g_na*I_iNa_m_na(n,:).^3.*I_iNa_h_na(n,:).*(I_V(n,:)-p.I_iNa_E_na))))+((((-p.I_iK_g_k*I_iK_m_k(n,:).^4.*(I_V(n,:)-p.I_iK_E_k))))+((((-((p.I_iPoisson_gext.*I_iPoisson_s(k,:))).*(I_V(n,:)-p.I_iPoisson_Eext))))+((((-p.I_E_iAMPA_gAMPA.*(I_E_iAMPA_s(n,:)*p.I_E_iAMPA_netcon).*(I_V(n,:)-p.I_E_iAMPA_EAMPA))))+((((-p.I_I_iGABAa_gGABAa.*(I_I_iGABAa_s(n,:)*p.I_I_iGABAa_netcon).*(I_V(n,:)-p.I_I_iGABAa_EGABAa))))))))));
  I_iE(n,:)=((((-((p.I_iPoisson_gext.*I_iPoisson_s(k,:))).*(I_V(n,:)-p.I_iPoisson_Eext))))+((((-p.I_E_iAMPA_gAMPA.*(I_E_iAMPA_s(n,:)*p.I_E_iAMPA_netcon).*(I_V(n,:)-p.I_E_iAMPA_EAMPA))))));
  I_iI(n,:)=((((-p.I_I_iGABAa_gGABAa.*(I_I_iGABAa_s(n,:)*p.I_I_iGABAa_netcon).*(I_V(n,:)-p.I_I_iGABAa_EGABAa)))));
  I_iL_i_l(n,:)=-p.I_iL_g_l*(I_V(n,:)-p.I_iL_E_l);
  I_iNa_i_Na(n,:)=-p.I_iNa_g_na*I_iNa_m_na(n,:).^3.*I_iNa_h_na(n,:).*(I_V(n,:)-p.I_iNa_E_na);
  I_iK_i_k(n,:)=-p.I_iK_g_k*I_iK_m_k(n,:).^4.*(I_V(n,:)-p.I_iK_E_k);

    n=n+1;
  end
end

T=T(1:downsample_factor:ntime);

end
