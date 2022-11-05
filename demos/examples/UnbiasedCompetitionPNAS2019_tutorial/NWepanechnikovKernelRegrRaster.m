function [inst_rate,t] = NWepanechnikovKernelRegrRaster(time,raster,pop,kwidth,Ts,flag_interp)

%  Nadaraya-Watson smoothing kernel regression of a rastergram using a Epanechnikov kernel

  if nargin < 6
    flag_interp = 0; % set this to 1 if you want to interpolate to the previous sampling frequency 1/(dsfact*dt)
  end

  if isempty(pop)
    pop = unique(raster(:,2));
  end

  cumSpikesPop_dt = zeros(size(time));
  for  iter_dt = 1:length(time)
    raster_dt = raster(raster(:,1)==time(iter_dt),:);
    cumSpikesPop_dt(iter_dt) = sum(ismember(raster_dt(:,2),pop));
  end
  interval = time(end)-time(1);
  ncells = length(pop);
  rate = sum(cumSpikesPop_dt)/interval/ncells;

  x = time;
  y = cumSpikesPop_dt;
  h = kwidth;
  % rcontrol = eksr(x,y,h);
  % rcontrol.fnorm = rcontrol.f/mean(rcontrol.f)*rate;

  % (down-)sampling spike times to Fs = 1/Ts (e.g. 1 kHz) because the kernel regression takes too long for small dt otherwise
  t = 0:Ts:max(time)+Ts;
  t = t(nearest(t,min(time))):Ts:t(nearest(t,max(time)));
  try
    dt = time(2)-time(1);
    cumSpikesPop_dt = resample(cumSpikesPop_dt,1,round(Ts/dt),2);
  catch
    time = double(time);
    dt = time(2)-time(1);
    cumSpikesPop_dt = resample(cumSpikesPop_dt,1,round(Ts/dt),2);
  end

%    if ~exist('t','var')
%      t = time;
%      flag_interp = 0;
%    end

  x = t;
  y = cumSpikesPop_dt;
  h = kwidth;
  r = eksr(x,y,h);
  r.fnorm = r.f/mean(r.f)*rate;
  % if min(r.fnorm) < 0
  %   r.fnorm = r.fnorm - min(r.fnorm);
  % end

  if ~exist('scalingFig','var')
    scalingFig = 2;
  end

  if flag_interp
    inst_rate = interp1(t,r.fnorm,time,'pchip');
    t = time;
  else
    inst_rate = r.fnorm;
    t = r.x;
  end

%    figure('color','none','visible','off')
%    hold on
%    set(gca,'layer','top','color','none')
%    plot(time,rcontrol.fnorm,'-','color',[0,0.6,1],'linewidth',1*scalingFig)
%    plot(t,r.fnorm,'-','color',[1,0.4,0.4],'linewidth',1*scalingFig)
%    title('Gaussian kernel regression','fontSize',16*scalingFig)
%    xlabel('Time (s)','fontSize',16*scalingFig)
%    ylabel('Instantaneous firing rate','fontSize',16*scalingFig)
%    set(gca,'fontSize',16*scalingFig,'LineWidth',1*scalingFig,'TickDir','out','Box','off')
%    plot2svg('gaussRegr.svg')
%
%    if flag_interp
%      figure('color','none','visible','off')
%      hold on
%      set(gca,'layer','top','color','none')
%      plot(time,inst_rate,'-','color',[0,0.6,1],'linewidth',1*scalingFig)
%      plot(t,r.fnorm,'-','color',[1,0.4,0.4],'linewidth',1*scalingFig)
%      title('Gaussian kernel regression vs interpolated curve','fontSize',16*scalingFig)
%      xlabel('Time (s)','fontSize',16*scalingFig)
%      ylabel('Instantaneous firing rate','fontSize',16*scalingFig)
%      set(gca,'fontSize',16*scalingFig,'LineWidth',1*scalingFig,'TickDir','out','Box','off')
%      plot2svg('gaussRegrVSinterpCurve.svg')
%
%      figure('color','none','visible','off')
%      hold on
%      set(gca,'layer','top','color','none')
%      plot(time,inst_rate-rcontrol.fnorm,'-','color',[1,0.4,0.4],'linewidth',1*scalingFig)
%      title('Diff Gaussian kernel regression','fontSize',16*scalingFig)
%      xlabel('Time (s)','fontSize',16*scalingFig)
%      ylabel('Instantaneous firing rate','fontSize',16*scalingFig)
%      set(gca,'fontSize',16*scalingFig,'LineWidth',1*scalingFig,'TickDir','out','Box','off')
%      plot2svg('diffGaussRegrVSinterpCurve.svg')
%    end
end
