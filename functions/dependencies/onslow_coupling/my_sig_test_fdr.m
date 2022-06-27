function [h,p] = my_sig_test_fdr(pacval, pac_shf_values, alpha, measure)
% function [h,p] = my_sig_test_fdr(pacval, pac_shf_values, alpha, measure)
%
% Returns a binary value for significance 'h' and a p-value 'p', given a
% value under analysis 'pacval', a distribution of values for comparison
% 'pac_shf_values' and a significance level 'alpha'
%
% INPUTS:
% pacval - PAC value 
% pac_shf_values - distribution of PAC values obtained from shuffled data
% sets 
% alpha - significance level, i.e. if alpha = 0.05 then 'pacval' must fall
% in the top 5% of the 'pac_shf_values' to be deemed significant
%
% Author: Angela Onslow, May 2010

if strcmp(measure, 'esc')
    
    p_tail1 = mean (pac_shf_values > pacval);
    p_tail2 = mean (pac_shf_values < pacval);
    
    if (p_tail1 <= alpha/2)
        h = 1;
        p = 2*(p_tail1);
    elseif (p_tail2 <= alpha/2)
        h = 1;
        p = 2*(p_tail2);
    elseif (p_tail1 > alpha/2)
        h = 0;
        p = 2*(p_tail1);
    elseif (p_tail2 > alpha/2)
        h = 0;
        p= 2*(p_tail2);
    end
 
    
else
    
    % percentage of times pacval is greater than a pac_shf_value
    p = mean (pac_shf_values > pacval);
    
    h = (p <= alpha);
    
end

if p == 0
    p =  (1/length(pac_shf_values))-(1/(10*length(pac_shf_values)));
end