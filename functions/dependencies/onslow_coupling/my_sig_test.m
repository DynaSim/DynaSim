function [h,p] = my_sig_test(pacval, pac_shf_values, alpha)
% function [h,p] = my_sig_test(pacval, pac_shf_values, alpha)
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
  

p = mean (pac_shf_values > pacval);

h = (p <= alpha);