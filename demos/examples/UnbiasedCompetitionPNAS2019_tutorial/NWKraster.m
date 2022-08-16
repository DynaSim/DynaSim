function [inst_rate,t] = NWKraster(time,raster,pop,kwidth,Ts,flag_interp,kernel)

%  Nadaraya-Watson smoothing kernel regression of a rastergram

  if nargin < 7
    kernel = 'E';
  end
  if nargin < 6
    flag_interp = 0; % set this to 1 if you want to interpolate to the previous sampling frequency 1/(dsfact*dt)
  end

  if isempty(pop)
    pop = unique(raster(:,2));
  end

  h = kwidth;
  ncells = numel(pop);

  dt = double(time(2)-time(1));
  if dt ~= Ts
    % (down-)sampling spike times to Fs = 1/Ts (e.g. 1 kHz) because the kernel regression takes too long for small dt otherwise
    t = min(time)-Ts:Ts:max(time)+Ts;
    t = t(nearest(t,min(time))):Ts:t(nearest(t,max(time)));
    cumSpikesPop_Ts = zeros(size(t));
    for  iter_Ts = 2:numel(t)
      raster_Ts = raster(raster(:,1) > t(iter_Ts-1) & raster(:,1) <= t(iter_Ts),:);
      cumSpikesPop_Ts(iter_Ts) = sum(ismember(raster_Ts(:,2), pop));
    end
    interval = t(end)-t(1);
    rate = sum(cumSpikesPop_Ts)/interval/ncells;
  else
    t = time;
    cumSpikesPop_Ts = cumSpikesPop_dt;
  end

  x = t;
  y = cumSpikesPop_Ts;

  if any(cumSpikesPop_Ts)
    switch kernel
      case 'E' % Epanechnikov
        r_Ts = eksr(x,y,h);
      case 'G' % Gaussian
        r_Ts = gksr(x,y,h);
      case 'L' % Laplacian
        r_Ts = lksr(x,y,h);
      otherwise
        disp('unknown kernel, using Epanechnikov')
        r_Ts = eksr(x,y,h);
    end

    if min(r_Ts.f) < 0
      warning('correcting negative rates')
      r_Ts.f(r_Ts.f < 0) = 0;
    end
    r_Ts.fnorm = r_Ts.f/mean(r_Ts.f)*rate;

    r = r_Ts;
    if flag_interp
      inst_rate = interp1(r.x,r.fnorm,time,'pchip');
      t = time;
    else
      inst_rate = r.fnorm;
      t = r.x;
    end
  else
    if flag_interp
      t = time;
    end
    inst_rate = zeros(size(t));
  end
end
