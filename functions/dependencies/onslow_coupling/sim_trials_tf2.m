function [trialmat trialmat2 avg_snr] = sim_trials_tf2(lf, hf, num_trials, Fs, num_sec, noise_ratio, amp_ratio)
% function [trialmat trialmat2 avg_snr] = sim_trials_tf2(lf, hf, num_trials, Fs, num_sec, noise_ratio, amp_ratio)
%
% This function creates 2 matrices of simulated trial data (trials stored as columns)
% which contain 'lf' modulating 'hf' PAC. Each trial is 'num_sec' seconds in length
% and the PAC occurs after a period of two seconds and lasts for a period of 1 to 2
% seconds; the rest of the signals' duration contains white noise. Each trial
% contains a random phase shift (<= pi/2) to distinguish it from the others.
% One matrix contains a brief period of PAC, the other contains a
% corresponding time period of the lower, modulating PAC frequency only.
%
% Inputs:
%
% lf: desired lower, modulating frequency
%
% hf: desired higher, modulated frequency
%
% num_trials: number of trials required 
%
% Fs: assumed sampling frequency for the trials
%
% num_sec: length of the desired signals in seconds
%
% noise_ratio: amplitude of white noise added to the PAC portion of the signal
%
% amp_ratio: amplitude ratio of the two components of the PAC signal
%
% Outputs:
%
% trialmat: matrix of simulated trials containing a brief period of PAC, 
% (Fs*num_sec) x num_trials in size
%
% trialmat2: matrix of simulated trials containing a brief period of lower,
% modulating PAC frequency, (Fs*num_sec) x num_trials in size 
% 
%
% Author: Angela Onslow, December 2009


trialmat = zeros((Fs*num_sec), num_trials);

lfvar = randi([-2 2], num_trials, 1);
hfvar = randi([-2 2], num_trials, 1);
phshift = pi/2.*rand(num_trials,1);

snr_sum =0;

for i = 1:num_trials
    
    
    [simsig simsigmod snr] = proc7((lf + lfvar(i)), (hf+hfvar(i)), Fs, num_sec, noise_ratio, amp_ratio, phshift(i));
    
    sigfirst = 0.5*randn((Fs*num_sec),1);
    sigsecond = 0.5*randn((Fs*num_sec),1);
    
    tshift = floor(Fs/2) + randi(floor(Fs/2), 1, 1);
    
    trialmat(:,i) = [sigfirst(1:(2*Fs));simsig(((2*Fs)+1):(tshift+1+(2*Fs)));sigfirst((tshift+2+(2*Fs)):(Fs*num_sec))];
    trialmat2(:,i) = [sigsecond(1:(2*Fs));simsigmod(((2*Fs)+1):(tshift+1+(2*Fs)));sigsecond((tshift+2+(2*Fs)):(Fs*num_sec))];
    
    snr_sum = snr_sum + snr;
    
end

avg_snr = snr_sum/num_trials;