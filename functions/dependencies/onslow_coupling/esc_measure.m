function escval = esc_measure(ph_sig, amp_sig, avg)
% function escval = esc_measure(ph_sig, amp_sig, avg)
%
% Returns a value (or vector of values) for the ESC measure calculated
% between two signals. Signals may contain one of more trials. Multiple
% trials may be averaged so as to return one ESC value or a vector of
% ESC values calculated for each trial may be returned, depending on the 
% 'avg' argument. Signals should be passed as column vectors, multiple 
% trials stored as multiple columns.
%
% INPUTS:
% ph_sig - signal filtered for a lower, modulating frequency (e.g. theta
% band oscillations)
%
% amp_sig - signal filtered for a higher, modulated frequency (e.g. gamma
% band oscillations)
%
% avg - string, either 'y' or 'n', determines whether ESC values are
% averaged over trials or returned as a vector
%
% Author: Angela Onslow, May 2010

escsum = 0;
if size(ph_sig, 2) ~= size(amp_sig, 2)
    sprintf('Error - Signals must have the same number of trials')
    return
end
num_trials = size(ph_sig, 2);

if strcmp(avg, 'y')
    
    %Average over trials using the Fisher transform
    for c = 1:num_trials
            
            r = corrcoef(ph_sig(:,c), amp_sig(:,c));
            escsum = escsum + atanh(r(1,2));
            
    end

    escsum = escsum/num_trials;
    escval = tanh(escsum);
    
else
    escval = zeros(num_trials,1);

    for i = 1:num_trials
        r = corrcoef(ph_sig(:,i), amp_sig(:,i));
        escval(i,1) = r(1,2);
    end
end
