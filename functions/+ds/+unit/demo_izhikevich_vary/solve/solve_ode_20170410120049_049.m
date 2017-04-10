function [T,pop1_v,pop1_u]=solve_ode
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
% ------------------------------------------------------------
% Initial conditions:
% ------------------------------------------------------------
% seed the random number generator
rng(p.random_seed);
t=0; k=1;

% STATE_VARIABLES:
pop1_v = zeros(nsamp,p.pop1_Npop);
  pop1_v(1,:) = -70*ones(1,p.pop1_Npop);
pop1_u = zeros(nsamp,p.pop1_Npop);
  pop1_u(1,:) = -20*ones(1,p.pop1_Npop);
% ###########################################################
% Numerical integration:
% ###########################################################
n=2;
for k=2:ntime
  t=T(k-1);
  pop1_v_k1=.04*pop1_v(n-1)^2+5*pop1_v(n-1)+140-pop1_u(n-1)+p.pop1_I;
  pop1_u_k1=p.pop1_a*(p.pop1_b*pop1_v(n-1)-pop1_u(n-1));
  t=t+.5*dt;
  pop1_v_k2=.04*pop1_v(n-1)^2+5*pop1_v(n-1)+140-pop1_u(n-1)+p.pop1_I;
  pop1_u_k2=p.pop1_a*(p.pop1_b*pop1_v(n-1)-pop1_u(n-1));
  pop1_v_k3=.04*pop1_v(n-1)^2+5*pop1_v(n-1)+140-pop1_u(n-1)+p.pop1_I;
  pop1_u_k3=p.pop1_a*(p.pop1_b*pop1_v(n-1)-pop1_u(n-1));
  t=t+.5*dt;
  pop1_v_k4=.04*pop1_v(n-1)^2+5*pop1_v(n-1)+140-pop1_u(n-1)+p.pop1_I;
  pop1_u_k4=p.pop1_a*(p.pop1_b*pop1_v(n-1)-pop1_u(n-1));
  % ------------------------------------------------------------
  % Update state variables:
  % ------------------------------------------------------------
  pop1_v(n)=pop1_v(n-1)+(dt/6)*(pop1_v_k1+2*(pop1_v_k2+pop1_v_k3)+pop1_v_k4);
  pop1_u(n)=pop1_u(n-1)+(dt/6)*(pop1_u_k1+2*(pop1_u_k2+pop1_u_k3)+pop1_u_k4);
  % ------------------------------------------------------------
  % Conditional actions:
  % ------------------------------------------------------------
  conditional_test=(pop1_v(n)>=30);
  if any(conditional_test), pop1_v(n,conditional_test)=p.pop1_c;pop1_u(n,conditional_test)=pop1_u(n,conditional_test)+p.pop1_d; end
  n=n+1;
end
T=T(1:downsample_factor:ntime);
