function [S,dynrate,spikeevents] = getPoissonDutyCycleGating(tOn,tOff,latency_,freq_,normFreqSigma_,phase_,widthSigma,baseline,dc_,ac_,slopeDutyCycle_,thresholdDutyCycle_,tau,kick,T,N)

  if nargin<1,  tOn = 0;                  end
  if nargin<2,  tOff = 0;                 end
  if nargin<3,  latency_ = 0;             end
  if nargin<4,  freq_ = 0;                end
  if nargin<5,  normFreqSigma_ = 0.03;    end
  if nargin<6,  phase_ = 0;               end
  if nargin<7,  widthSigma = 0.001;       end % 0.001 represents an abrupt connectivity transition (in contrast to 0.1; it only applies to dc+ac not the baseline)
  if nargin<8,  baseline = 0;             end % Hz
  if nargin<9,  dc_ = 0;                  end % Hz
  if nargin<10, ac_ = 0;                  end % Hz
  if nargin<11, slopeDutyCycle_ = 1;      end % in ms
  if nargin<12, thresholdDutyCycle_ = 0;  end
  if nargin<13, tau = 2;                  end
  if nargin<14, kick = 1;                 end
  if nargin<15, T = (0:.01:1000)';        end
  if nargin<16, N = 1;                    end

  dt=0.5*(T(2)-T(1)); % RK steps of 0.5dt
  interval=T(end)-T(1);

  numStim = length(tOn(~isnan(tOn)));

  % broadcasting sizes with respect to number of stimuli
  if numStim > size(freq_,1), freq = 1e-3*ones(numStim,1)*freq_; else, freq = 1e-3*freq_; end % in kHz
  if numStim > size(normFreqSigma_,1), normFreqSigma = ones(numStim,1)*normFreqSigma_; else, normFreqSigma = normFreqSigma_; end
  if numStim > size(phase_,1), phase = ones(numStim,1)*phase_; else, phase = phase_; end
  if numStim > size(dc_,1), dc = ones(numStim,1)*dc_; else, dc = dc_; end
  if numStim > size(ac_,1), ac = ones(numStim,1)*ac_; else, ac = ac_; end
  if numStim > size(slopeDutyCycle_,1), slopeDutyCycle = ones(numStim,1)*slopeDutyCycle_; else, slopeDutyCycle = slopeDutyCycle_; end
  if numStim > size(thresholdDutyCycle_,1), thresholdDutyCycle = ones(numStim,1)*thresholdDutyCycle_; else, thresholdDutyCycle = thresholdDutyCycle_; end
  if numStim > size(latency_,1), latency = ones(numStim,1)*latency_; else, latency = latency_; end

  time = 0:dt:interval;
  S = zeros(N,length(time));
  S_ini = zeros(N,1);
  dynrate = baseline*ones(N,length(time));
  ratecomp = zeros(numStim,length(time));
  for i = 1:numStim
    timeWindow = zeros(1,length(time));
    if latency(i) ~= 0
      timeWindow(time >= tOn(i)) = 1;
    else
      timeWindow(time >= tOn(i) & time <= tOff(i)) = 1;
    end
    if freq(i)
      % freq modulation
      numFreqs = 1000;
      step = 2*pi/numFreqs;
      Ph = load('sinPhases'); % from test.m in this directory
      % phase_ic = 0:step:2*pi*(1-1/numFreqs);
      % shuffle = randperm(numFreqs);
      % Ph.phase_ic = phase_ic(shuffle);
      % save('sinPhases','Ph')
      freqSigma = normFreqSigma(i)*freq(i);
      freqSet = -freq(i)/5:2*freq(i)/5/(numFreqs-1):freq(i)/5;
      freqVar = exp(-freqSet.^2/(2*freqSigma^2));
      m = sum(freqVar)/numFreqs;
      sumSin = zeros(size(time));
      for j = 1:numFreqs
        sumSin = sumSin + freqVar(j)*cos(2*pi*freqSet(j)*time + Ph.phase_ic(j));
      end

      sinSignal = sin(2*pi*freq(i)*time-m*sumSin+phase(i));
      % sinSignal = sin(2*pi*freq(i)*time-pi/2+phase(i));

      dutySquareSignal = sign(sinSignal-thresholdDutyCycle(i));

      ups_tmp_4coder = 1+find(diff(dutySquareSignal)>0);
      ups_tmp2_4coder = 1+find(diff(ups_tmp_4coder)==1);
      % ups_tmp_4coder(ups_tmp2_4coder) = []; %  coder does not seem to like this
      ups_tmp_4coder(ups_tmp2_4coder) = 0;
      ups_tmp3_4coder = ups_tmp_4coder(ups_tmp_4coder>0);
      ups_tmp3_4coder = time(ups_tmp3_4coder)';
      if ~isempty(ups_tmp3_4coder)
        ups = [ups_tmp3_4coder(1)-1/freq(i);ups_tmp3_4coder;ups_tmp3_4coder(end)+1/freq(i)];
      else
        ups = [];
      end

      downs_tmp_4coder = 1+find(diff(dutySquareSignal)<0);
      downs_tmp2_4coder = 1+find(diff(downs_tmp_4coder)==1);
      % downs_tmp_4coder(downs_tmp2_4coder) = []; %  coder does not seem to like this
      downs_tmp_4coder(downs_tmp2_4coder) = 0;
      downs_tmp3_4coder = downs_tmp_4coder(downs_tmp_4coder>0);
      downs_tmp3_4coder = time(downs_tmp3_4coder)';
      if ~isempty(downs_tmp3_4coder)
        downs = [downs_tmp3_4coder(1)-1/freq(i);downs_tmp3_4coder;downs_tmp3_4coder(end)+1/freq(i)];
      else
        downs = [];
      end

      perDutySigmoidSignal = sum(1./(1+exp(-(ones(size(ups))*time-ups*ones(size(time)))/slopeDutyCycle(i))),1) - sum(1./(1+exp(-(ones(size(downs))*time-downs*ones(size(time)))/slopeDutyCycle(i))),1);
      normPerDutySigmoidSignal = (perDutySigmoidSignal - min(perDutySigmoidSignal))/(max(perDutySigmoidSignal) - min(perDutySigmoidSignal));
      ratecomp(i,:) = dc(i) + ac(i)*2*(normPerDutySigmoidSignal - 0.5);
      ratecomp(i,ratecomp(i,:)<0) = 0;
    else
      ratecomp(i,:) = dc(i);
    end
    ratecomp(i,:) = ratecomp(i,:).*timeWindow;
    if latency(i) ~= 0
      ratecomp(i,time>=tOn(i)) = ratecomp(i,time>=tOn(i)).*(1-exp(-(time(time>=tOn(i))-tOn(i))/latency(i)));
      ratecomp(i,time>=tOff(i)) = ratecomp(i,time>=tOff(i)).*exp(-(time(time>=tOff(i))-tOff(i))/latency(i));
    end
  end
  % x_N = ((1:N)-0.5)/N;
  % x_c = 0.25; % center of first category
  % conn = 1./(1+exp(-cos(2*pi*(x_N-x_c))/widthSigma));
  % dynrate = dynrate + conn'*sum(ratecomp,1);
  dynrate = dynrate + ones(N,1)*sum(ratecomp,1);
  % figure
  % plot(time,dynrate(1,:))
  if any(dynrate(1,:) < 0)
    error('dynrate is < 0')
  end
  % if any(freq) > 0
  %   plotVariable(time,dynrate(1,:),[],[],[],0,0,unique([tOn,tOff]),1,'Time(ms)','Rate (Hz)','k',[],'Plot dynrate','../../dynrate_m');
  %   error('stop')
  % end
  % plot(time, dynrate(1,:))
  % hold on
  % drawnow
  [s,spikeevents] = nonhomPoissonGatingGenerator(S_ini,dynrate/1e3,tau,kick,N,interval,dt); % dynrate in kHz
  S=s';
end
