function shf_sig = shuffle_trials(sig)
% function shf_sig = shuffle_trials(sig)
%
% Shuffles the trials the columns of sig (assumed to be trials) randomly
%
% INPUTS:
%
% sig - a matrix or cell array of matrices
%
% OUTPUTS:
% shf_sig - either a vector or a cell array depending on the input 'sig'.
% When used within find_pac_shf.m this output is of the same dimension as 
% filt_sig_mod and filt_sig_pac: number of cells - number of frequency bins
% and each cell element is a matrix(num_data_points, num_trials)
%
% Author: David Stanley, Boston University, March 2015


    sig_type = class(sig);

    switch sig_type
        case 'double'
            shf_sig = shuffle_trials_matrix(sig);
        case 'cell'
            shf_sig = cellfun(@shuffle_trials_matrix,sig,'UniformOutput',0);
    end

end


function shf_sig = shuffle_trials_matrix(sig)
    
    Ntrials = size(sig,2);
    ind = randperm(Ntrials);
    shf_sig = sig(:,ind);
end
                