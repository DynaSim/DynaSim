function [raster,rndRaster]  = computeRaster(t,V,threshold)
  if nargin<3, threshold=0; end
  dt = t(2)-t(1);
  raster = [];
  [indTimes,neuronSpikes] = find (V > threshold);
  if ~isempty(neuronSpikes)
    tSpikes = t(indTimes); % in s
    raster(:,1) = tSpikes;
    raster(:,2) = neuronSpikes;
    raster(diff(raster(:,2))==0 & diff(raster(:,1))<=1.05*dt,:) = []; % removing artificial spikes that come from two consecutive voltages above 0 mV
    [~,indSort] = sort(raster(:,1));
    raster = raster(indSort,:);
  end

  if ~isempty(raster)
    rndInd = randperm(max(raster(:,2)));
    rndRaster = raster;
    rndRaster(:,2) = rndInd(rndRaster(:,2));
  else
    rndRaster = raster;
  end
