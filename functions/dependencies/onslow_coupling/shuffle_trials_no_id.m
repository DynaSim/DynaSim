function shf_sig = shuffle_trials_no_id(sig)
% function shf_sig = shuffle_trials2(sig)
%
% Shuffles the trials the columns of sig (assumed to be trials) randomly,
% ensuring that none of the trials remain in their original positions. For
% example, this function could map [1 2 3 4] onto [2 1 4 3] but not onto
% [2 1 3 4]. In contrast, shuffle_trials could potentially do both.
%
% INPUTS:
%
% sig - Input signal in the form of a matrix or cell array of matrices.
% Columns are trials and rows are samples.
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
            shf_sig = shuffle_trials_matrix_for(sig);
        case 'cell'
            shf_sig = cellfun(@shuffle_trials_matrix_for,sig,'UniformOutput',0);
    end

end



function shf_sig = shuffle_trials_matrix_for(sig)           % Faster version without identity, ironically, using a for loop
    % shuffle_trials_matrix_for
    % Get all possible pairs
    
    Ntrials = size(sig,2);
    ind = randperm_no_identity(Ntrials);
    shf_sig = sig(:,ind);
end


function ind = randperm_no_identity(N)
    % randperm_no_identity
    % Works just like matlab's randperm(N) command, returning a random
    % permutation of the values 1:Ntrials. However, it will not return any
    % identity values.
    % For example, randperm(4) could return [2 1 4 3] or [2 1 3 4]
    % But randperm_no_identity(4) could only return [2 1 4 3] because the
    % array [2 1 3 4] contains an identity (the 3,4 part).
    % Uses a montecarlo type approach on the original matlab randperm
    % command
    % To do : modify to work with randperm(N,k)
%     tic
    ind = randperm(N);
    ind_id = 1:N;
    bads = find((ind - ind_id) == 0);       % Find the indicies of identities
    
    while ~isempty(bads)    % If there are any identities in the random permutation
        shuff_ind = floor(unifrnd(1,N+1-0.1,1,1));          % Find some other index at random
        ind([bads(1) shuff_ind]) = ind([shuff_ind bads(1)]);    % Swap the bad index with the other index, and hope this removes the identity; if not, we will repeat
        bads = find((ind - ind_id) == 0);           % Recalculate the number of identities remaining.
    end
%     ind;
%     toc
    

end


% OLD CODE AND FUNCTIONS

function shf_sig = shuffle_trials_matrix_nchoosek(sig)           % Shuffle trials without identity . This is very slow for huge number of trials
    % shuffle_trials_matrix_nchoosek
    % Get all possible pairs
    
    Ntrials = size(sig,2);
    tic
    %Ntrials = 900;
    pairs = nchoosek(1:Ntrials,2);                  
    pairs = [pairs; fliplr(pairs)];                 % Double for symmetry
    pairs = sortrows(pairs);
    
    % Now we want to pick 1 pair for each of the Ntrials phase vectors
    pairs = reshape(pairs',[2,Ntrials-1,Ntrials]);
    
    i=floor(unifrnd(1,Ntrials-0.01,1,Ntrials));     % We want this to be from 1 to Ntrials-1. Subtract 0.001 and floor takes care of the rest.
    j=1:Ntrials;
    
    sz = size(pairs); sz=sz(2:3);
    ind = sub2ind(sz,i,j);
    
    pairs_chosen = pairs(:,ind);
    pairs_chosen = pairs_chosen';
    toc                         % Takes 2.3 seconds for 900 trials
    
    shf_sig = sig(:,pairs_chosen(:,2));
end
         


function ind = randperm_no_identity_broken(Ntrials)        % Fast but doesn't work - can get stuck on the last digit
    % randperm_no_identity_broken
    
    tic
    ind = zeros(1,Ntrials);
    pool = 1:Ntrials;
    for i = 1:Ntrials
        
        chosen = floor(unifrnd(1,length(pool)+1-0.1,1,1));  % Choose a random value from the pool
        while pool(chosen) == i
            chosen = floor(unifrnd(1,length(pool)+1-0.1,1,1));  % If it's an identity, rechoose
        end
        
        ind(i) = pool(chosen);
        pool = pool([1:chosen-1, chosen+1:end]);
    end
    ind;
    
    toc

end

