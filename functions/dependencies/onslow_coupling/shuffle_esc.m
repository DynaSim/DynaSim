function shf_sig = shuffle_esc(sig, Fs)
% function shf_sig = shuffle_esc(sig, Fs)
%
% This function shuffles a signal using random insertion type shuffling. 
% The signal 'sig' is divided into a number of sections, equal to the
% number of seconds of data or 1000, which ever is smaller. These sections 
% are of random lengths, chosen from a uniform distribution. The sections 
% are then randomly re-ordered. 
% The input signal 'sig' must contain several seconds on data (>2) and must
% be a vector (i.e. number of trials = 1) or a cell array containing 
% multiple signals
%
% INPUTS:
% sig - input signal which is to be shuffled, passed as a column vector
% Fs - sampling frequency of 'sig', in Hz
%
% OUTPUTS:
% shf_sig - either a vector or a cell array depending on the input 'sig'.
% When used within find_pac_shf.m this output is of the same dimension as 
% filt_sig_mod and filt_sig_pac: number of cells - number of frequency bins
% and each cell element is a matrix(num_data_points, num_trials)
%
% Author: Rafal Bogacz, Angela Onslow, May 2010

sig_type = class(sig);

if iscell(sig)
    ybins = size(sig,1);
    xbins = size(sig,2);
    shf_sig_cell = cell(size(sig));
    num_sec = ceil((size(sig{1,1},1))/Fs);
else
    shf_sig = [];
    num_sec = ceil((size(sig,1))/Fs);
end

if num_sec > 1000
    num_sec = 1000;
end


switch sig_type
    
    case 'double'
        
        % Choose num_sec random 'cut' positions
        dpsplit = ceil(size(sig,1).*rand(num_sec,1));
        % Arrange these in ascending order
        dpsplit = sort (dpsplit);

        start(1)=1;
        start(2:num_sec)=dpsplit(1:num_sec-1);
        ending(1:num_sec-1)=dpsplit(1:num_sec-1)-1;
        ending(num_sec) =  size(sig,1);

        order = randperm(num_sec);
    
        for c = 1:num_sec
        
            %shuffle the signal
            shf_sig = [shf_sig; sig(start(order(c)):ending(order(c)),:)];
        
        end


    case 'cell'
        for i = 1:ybins
            for j = 1:xbins
                
                current_sig = sig{i,j};
                
                % Choose num_sec random 'cut' positions
                dpsplit = ceil(size(sig{1,1},1).*rand(num_sec,1));
                % Arrange these in ascending order
                dpsplit = sort (dpsplit);

                start(1)=1;
                start(2:num_sec)=dpsplit(1:num_sec-1);
                ending(1:num_sec-1)=dpsplit(1:num_sec-1)-1;
                ending(num_sec) =  size(sig{1,1},1);

                order = randperm(num_sec);
    
                for c = 1:num_sec
        
                    %shuffle the signal
                    shf_sig_cell{i,j} = [shf_sig_cell{i,j}; current_sig(start(order(c)):ending(order(c)),:)];
        
                end
                
            end
        end
        
        shf_sig = shf_sig_cell;
end
                