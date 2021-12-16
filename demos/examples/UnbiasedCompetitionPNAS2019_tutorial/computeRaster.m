function [raster,rndRaster]  = computeRaster(t,V)
  dt = t(2)-t(1);
  raster = [];
  [indTimes,neuronSpikes] = find (V > 0);
  if ~isempty(neuronSpikes)
    tSpikes = t(indTimes); % in s
    raster(:,1) = tSpikes;
    raster(:,2) = neuronSpikes;
    [~,indSortN] = sort(raster(:,2));
    raster = raster(indSortN,:);
    raster(diff(raster(:,2)) < 0.5 & diff(raster(:,1)) < 1.5*dt,:) = []; % removing artificial spikes that come from two consecutive voltages above 0 mV
    [~,indSortT] = sort(raster(:,1));
    raster = raster(indSortT,:);
  end

  if ~isempty(raster)
    rndInd = randperm(max(raster(:,2)));
    rndRaster = raster;
    rndRaster(:,2) = rndInd(rndRaster(:,2));
  else
    rndRaster = raster;
  end
end
