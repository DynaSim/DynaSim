function [trialmat trialmatmod] = sim_trials(lf, hf, num_trials, Fs, num_sec)
%function trialmat = sim_trials(lf, hf, num_trials, Fs, num_sec)
%
%This function creates a matrix of simulated trial data (trials stored as columns)
%which contain 'lf' modulating 'hf' PAC. Each trial is 'num_sec' seconds in length
%and the PAC occurs throughout the entire length of signal. Each trial
%contains a random phase shift (<= pi/2) to distinguish it from the others
%
%Inputs:
%
%lf: desired lower, modulating frequency
%
%hf: desired higher, modulated frequency
%
%num_trials: number of trials required 
%
%Fs: assumed sampling frequency for the trials
%
%num_sec: length of the desired signals in seconds
%
%Outputs:
%
%trialmat: matrix of simulated trials, (Fs*num_sec) x num_trials in size
%
%Example usage:
%
%trialmat = sim_trials(4, 60, 16, 1017, 10)
%
%Author: Angela Onslow, December 2009

trialmat = zeros((Fs*num_sec), num_trials);

lfvar = randi([-2 2], num_trials, 1);
hfvar = randi([-2 2], num_trials, 1);
phshift = pi/2.*rand(num_trials,1);


for i = 1:num_trials
    
    
    [trialmat(:,i) trialmatmod(:,i)] = proc7((lf + lfvar(i)), (hf+hfvar(i)), Fs, num_sec, 0.1, 1, phshift(i));
       
    
end