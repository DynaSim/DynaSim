function rms_out = rms_val(sig)
% function rms_out = rms_val(sig)
%
% Returns the root mean square value of the signal passed as 'sig'
%
% Author: Angela Onslow, May 2010

rms_out = norm(sig)/sqrt((length(sig)));

end