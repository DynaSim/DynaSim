function y = dlPowerSpectrumRatio(dlObj, opts)

% opts : struct with fields -> <lf1,hf1,lf2,hf2>

    dlPotentialIndices = contains(dlObj.dlVariables, '_V');
    dlPotentialIndices(1) = 1;
    dlPotentials = dlObj.dlOutputs(dlPotentialIndices);
    dlLabels = dlObj.dlVariables(dlPotentialIndices);

%     t = dlPotentials{1, 1};
    n = size(dlPotentials, 2);
    y = zeros(n, 1); 
    m = ceil(n/6);
            
    dtf = ceil(1 / (dlObj.dldT*dlObj.dlDownSampleFactor));
    
    lf1 = opts.lf1*dtf;
    hf1 = opts.hf1*dtf;
    lf2 = opts.lf2*dtf;
    hf2 = opts.hf2*dtf;
    
    for k = 1:m
        
%         freqCap = 0;

        for i = (k-1)*6+1:min((k*6), 6)

            x = dlPotentials{1, i+1};
%             fqs = linspace(1, 500, max(size(x)));
%             subplot((min(k*6, n-1) - (k-1)*6), 1, mod(i-1, (min(k*6, n-1) - (k-1)*6))+1);
            ffts = abs(fft(mean(x, 2))) * min(size(x)) / 1000;
            
            yf1 = smooth(ffts(lf1:hf1));
            yf2 = smooth(ffts(lf2:hf2));
            y(i) = mean(yf1) / mean(yf2);
%             fprintf("---> %s 's power ratio [%.2f-%.2f]Hz/[%.2f-%.2f]Hz = %.3f\n", dlLabels{1, i+1}, opts.lf1, opts.hf1, opts.lf2, opts.hf2, y(i));
            
%             area(fqs(lf1:hf1), yf1);grid("on");
%             area(fqs(lf2:hf2), yf2);grid("on");

%             if freqCap == 0
%                 freqCap = max(ffts(lf1:hf1))*1.2;
%                 ylim([0, freqCap]);
%             else
%                 ylim([0, freqCap]);
%             end
% 
%             ylabel(dlLabels(i+1));

        end
    end
    
%     fprintf("\n->Power spectrum ratios.");

end

