function [pacmat,freqvec_ph_final2,freqvec_amp_final] = cfc_measure(sig1,...
    sig2, avg, freqvec_ph, freqvec_amp, Fs, nfft, width, waitbar)
% 
% This function calculates the CFC measure when passed two unfiltered
% signals.
% 
% INPUTS:
% sig1 - signal containing (or believed to contain) PAC
% sig2 - signal containing (or believed to contain) the modulating signal
% avg - (string) determines whether CFC values should be averaged over, if
% == 'y' then CFC values will be averaged over trials
% freqvec_ph - vector of lower, modulating frequencies aimed to filter at
% (approximate due to the use of the mscohere function)
% freqvec_amp - vector of higher, modulated frequencies to filter at 
% (achieved via convolution with complex Morlet wavelets, see ENERGYVEC)
% Fs - sampling frequency, in Hz
% nfft - the number of points in FFT (applied by the mscohere function)
% width - width of the wavelet filter (defines the number of cycles of the 
% mother wavelet, see MORLET_WAV)
% waitbar - (string) determines whether to output a progress level to the
% command line, if == 'y' then a percentage value indicating the time
% remaining is displayed
%
% OUTPUTS:
% pacmat - matrix of PAC values (specifically CFC values)
% freqvec_ph_final2 - vector of centre frequencies which 'sig2' has been
% filtered at (x-axis of the PACgram)
% freqvec_amp_final - vector of centre frequencies which 'sig1' has been
% filtered at (y-axis of the PACgram)
%
% Author: Angela Onslow, May 2010



nfft = Fs;
% if mod(nfft,2) ~= 0
%     nfft = nfft+1;
% end


xbins = ceil((max(freqvec_ph) - min(freqvec_ph))/(diff(freqvec_ph(1:2))));
cent_freq_vec = zeros(xbins,1);

for i =1:xbins
    upper_bin = min(freqvec_ph)+i*(diff(freqvec_ph(1:2)));
    lower_bin = upper_bin-(diff(freqvec_ph(1:2)));
    cent_freq_vec(i) =  lower_bin + floor((upper_bin- lower_bin)/2);
end




freqvec_amp_final = freqvec_amp';

if mod(nfft,2) == 0
    length_Cxy = (nfft/2 +1);
else
    length_Cxy = (nfft+1)/2;
end

    
CRF = zeros(length(freqvec_amp), length_Cxy);
CRFcell = cell(length(freqvec_amp), length_Cxy);

if size(sig1, 2) ~= size(sig2, 2)
    sprintf('Error - Signals must have the same number of trials')
    return
end
num_trials = size(sig1, 2);




for c = 1:num_trials
    
    for k=1:length(freqvec_amp)
        E1 = energyvec(freqvec_amp(k),sig1(:,c)',Fs,width);
        if length(sig1) > (2*Fs)-1
            [CRF(k,:),freqvec_ph_final] = mscohere(E1',sig2(:,c),hanning(nfft),[],nfft,Fs);
        else
            [CRF(k,:),freqvec_ph_final] = mscohere(E1',sig2(:,c),[],[],nfft,Fs); 
        end
       
        if waitbar ==1
            if k ==1
                fprintf('%03i%% ', floor((k/length(freqvec_amp))*100));
            else
                fprintf('\b\b\b\b\b%03i%% ', floor((k/length(freqvec_amp))*100));
            end
            if k == length(freqvec_amp)
                fprintf('\n');
            end
        end
    end
    % N.B - divisions of freqVec1 = ceil((Fs/2)/((nfft/2) + 1));
    % freqVec1 will always contain (nfft/2) +1 elements if nfft is even
    % or (nfft+1)/2 elements if nfft is odd
    
    for i = 1:length(freqvec_amp)
        for j = 1:length_Cxy
            
            CRFcell{i,j} = [CRFcell{i,j} CRF(i,j)];
        end
    end
    
    
end


if strcmp(avg, 'y')
    CRF = cellfun(@mean, CRFcell);
      
    
    CRF = CRF(:,1:floor(end/2));
    freqvec_ph_final = freqvec_ph_final(1:ceil(end/2));
    
    % Get rid of first element correspoding to zero frequency and extra
    % elements which are calculated due to the mscohere method
    
    [C,ia] = intersect(freqvec_ph_final,cent_freq_vec);
       
    for c = 1:length(ia)
        freqvec_ph_final2(c) = freqvec_ph_final(ia(c));
        CRF2(:,c) = CRF(:,ia(c));
    end
    
else
    
    
    CRF = CRFcell;
    CRF = CRF(:,1:floor(end/2));
    freqvec_ph_final = freqvec_ph_final(1:floor(end/2));
    
    [C,ia] = intersect(freqvec_ph_final,cent_freq_vec);
       
    for c = 1:length(ia)
        freqvec_ph_final2(c) = freqvec_ph_final(ia(c));
        CRF2(:,c) = CRF(:,ia(c));
    end

end


pacmat = CRF2;
