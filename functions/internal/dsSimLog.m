function tstart = dsSimLog(k,logTimes,T,tstart)
%dsSimLog - if called, solve_ode files will inform about partial elapsed time at regular periods of the simulation
%
% Author: Salva Ardid
% Copyright (C) 2021 Salva Ardid

if tstart == 0
  tstart = tic;
elseif any(k == logTimes)
  elapsedTime = toc(tstart);
  elapsedTimeMinutes = floor(elapsedTime/60);
  elapsedTimeSeconds = rem(elapsedTime,60);
  if elapsedTimeMinutes
    logMS = 'Processed %g of %g ms (elapsed time: %g m %.3f s)\n';
    fprintf(logMS,T(k),T(end),elapsedTimeMinutes,elapsedTimeSeconds);
  else
    logS = 'Processed %g of %g ms (elapsed time: %.3f s)\n';
    fprintf(logS,T(k),T(end),elapsedTimeSeconds);
  end
end
