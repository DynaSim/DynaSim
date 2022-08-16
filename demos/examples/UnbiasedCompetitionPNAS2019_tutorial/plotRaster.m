function rate = plotRaster(rst,pool,tl,display_flag)

  if nargin < 3 || (isempty(pool) | isempty(tl))
    error('you should supply at least 3 input arguments, a cell array for the raster/-s, a dimensionally consistent cell array pool, and the time limits: [ti, tf]. Optionally, you can enable/disable the display of the raster plot. Default is display on. If you are only interested in the avg. FR computation, set the fourth input argument to false/off/0.')
  end

  rate = 0;
  if isempty(rst)
    return;
  end

  if nargin < 4
    display_flag = true;
  end

  if display_flag
    lineWidth = 1;
    fontSize = 16;
    colors = [
               0.0000    0.4470    0.7410
               0.8500    0.3250    0.0980
               0.9290    0.6940    0.1250
               0.4940    0.1840    0.5560
               0.4660    0.6740    0.1880
               0.3010    0.7450    0.9330
               0.6350    0.0780    0.1840
               0.0000    0.0000    0.0000
               33/255    113/255   181/255
               239/255   59/255    44/255
             ];
    figure % ('visible','on') %% this interferes with the live script spirit
    hold on
    set(gca,'layer','top')
  end

  npopulations = numel(rst);
  for i = 1:npopulations
    Npool(i) = numel(pool{i});
  end

  for ipopulation = 1:npopulations
    raster = rst{ipopulation};
    if ~isempty(raster)
      tSpikes = raster(:,1);
      raster = raster(tSpikes >= tl(1) & tSpikes <= tl(2),:);
      neuronTag = raster(:,2);
      rate(ipopulation) = sum(raster(ismember(neuronTag,pool{ipopulation}),1) >= tl(1) & raster(ismember(neuronTag,pool{ipopulation}),1) <= tl(2))/(1e-3*diff(tl))/Npool(ipopulation); % in sp/s
    else
      rate(ipopulation) = 0;
    end

    if display_flag
      plot(raster(:,1),raster(:,2),'ko','markerfacecolor',colors(ipopulation,:),'markerSize',4,'linewidth',lineWidth)
      if max(Npool) >= 75
        set(gca,'fontSize',fontSize,'LineWidth',lineWidth,'TickDir','out','Box','on','YTick',0:25:max(Npool));
      elseif max(Npool) >= 30
        set(gca,'fontSize',fontSize,'LineWidth',lineWidth,'TickDir','out','Box','on','YTick',0:10:max(Npool));
      else
        set(gca,'fontSize',fontSize,'LineWidth',lineWidth,'TickDir','out','Box','on','YTick',0:max(Npool)/2:max(Npool));
      end
        if ~exist('rate_label','var')
        rate_label = ['FR_', num2str(ipopulation), ' = ', num2str(rate(ipopulation),2),' sp/s'];
        else
        rate_label = [rate_label, ', FR_', num2str(ipopulation), ' = ', num2str(rate(ipopulation),2),' sp/s'];
        end
    end
  end

  if display_flag
    xlim(tl)
    if max(Npool) > 1
      ylim([0,max(Npool)+1])
    else
     ylim([0.5 1.5])
    end
    xlabel('Time (ms)','fontSize',fontSize)
    ylabel('Neurons','fontSize',fontSize)
    if ~isempty(rate_label)
        title(rate_label,'fontSize',fontSize,'FontWeight','Normal')
    end
  end
end
