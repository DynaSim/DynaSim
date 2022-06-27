function trialmat = sim_trials_tf(reg_f, diff_f, num_trials, Fs, num_sec, noise_ratio, amp_ratio1, amp_ratio2)
% function trialmat = sim_trials_tf(reg_f, diff_f, num_trials, Fs, num_sec, noise_ratio, amp_ratio1, amp_ratio2)
%
% This function creates a single matrix of trial data, with each trial as a
% separate column. Each trial containes power at 'reg_f' Hz + some noise, 
% after 2 seconds this main frequency component is changed to 'diff_f' Hz
% for a period of between 1 and 2 seconds (slightly different for each 
% trial) and then reverts back to 'reg_f' Hz again.
%
% Inputs:
%
% reg_f: desired regular frequency component, in Hz
%
% diff_f: desired different frequency component, in Hz 
%
% num_trials: number of trials required 
%
% Fs: assumed sampling frequency for the trials
%
% num_sec: length of the desired signals in seconds
%
% noise_ratio: amplitude of white noise to add to the signal
%
% amp_ratio1: amplitude of the regular portion of signal
%
% amp_ratio2: amplitude of the different portion of signal
%
% Outputs:
%
% trialmat: matrix of simulated trials, (Fs*num_sec) x num_trials in size
%
% Example usage:
%
% trialmat = sim_trials(10,4,3,1000,10,0.1,1,1);
%
% Author: Angela Onslow, May 2011

if nargin < 7
    amp_ratio1 = 1;
    amp_ratio2 = 1;
end

trialmat = zeros((Fs*num_sec), num_trials);

phshift = pi/2.*rand(num_trials,1);
phshift2 = pi/2.*rand(num_trials,1);

for i = 1:num_trials
    
    noise_sig1 = noise_ratio.*randn((Fs*num_sec),1);
    noise_sig2 = noise_ratio.*randn((Fs*num_sec),1);
    
    simsig = ((sin((reg_f/Fs)*2*pi*[1:(Fs*num_sec)]+phshift(i))).*amp_ratio1)';
    simsig2 = ((sin((diff_f/Fs)*2*pi*[1:(Fs*num_sec)]+phshift2(i))).*amp_ratio2)';
    
    simsig = simsig+noise_sig1;
    simsig2 = simsig2+noise_sig2;
    

    tshift = floor(Fs/2) + randi(floor(Fs/2), 1, 1);
    
    trialmat(:,i) = [simsig(1:(2*Fs));simsig2(((2*Fs)+1):(tshift+1+(2*Fs)));simsig((tshift+2+(2*Fs)):(Fs*num_sec))];
    
end