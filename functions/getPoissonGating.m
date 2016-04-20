function sext=getPoissonGating(baseline,dc,ac,freq,phase,onset,offset,tau,T,N,kernel,kick,ramp_dc_flag,ramp_ac_flag)

% default parameters
if nargin<1, baseline=100; end % Hz
if nargin<2, dc=1000; end % Hz
if nargin<3, ac=0; end % Hz
if nargin<4, freq=0; end % Hz
if nargin<5, phase=0; end % radians
if nargin<6, onset=T(1); end % ms
if nargin<7, offset=T(end); end % ms
if nargin<8, tau=2; end % ms
if nargin<9, T=(0:.01:1000)'; end % ms
if nargin<10, N=1; end % number of target cells
if nargin<11, kernel=ones(1,N); end % connectivity to target cells
if nargin<12, kick=1; end
if nargin<13, ramp_dc_flag=0; end % whether to ramp dc from 0 to dc over [onset,offset], else step at onset
if nargin<13, ramp_ac_flag=0; end % whether to ramp ac from 0 to ac over [onset,offset], else step at onset

dt=T(2)-T(1);
interval=T(end)-T(1);
latency=.1; 
fspread=.03; 
consigma=.001; % not used
s=getGenExtPoissonTotalGating(onset,offset,latency,freq/1000,fspread,phase,consigma,baseline/1000,dc/1000,ac/1000,tau,kick,N,interval,dt,kernel',ramp_dc_flag,ramp_ac_flag);
sext=s';
