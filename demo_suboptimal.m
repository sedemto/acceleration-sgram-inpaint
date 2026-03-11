clc;clear;close all

ltfatstart;
%% load signal
[signal,fs] = audioread("a08_violin.wav");

%% testing suboptimally shortened spectrograms

% most common case of overlap
w = 2048;
SNR_case1 = getSNRsuboptimal(signal,w,w/4,w,'hann','timeinv');

% another case
w2 = 2049;
SNR_case2 = getSNRsuboptimal(signal,w2,w2/3,w2,'hann','timeinv');

% another case
SNR_case3 = getSNRsuboptimal(signal,w,w/2,w,'hann','timeinv');

disp("SNR when w=2048, a=512: "+SNR_case1);
disp("SNR when w=2049, a=683: "+SNR_case2);
disp("SNR when w=2048, a=1024: "+SNR_case3);



function [SNR] = getSNRsuboptimal(signal,w,a,M,wtype,phasetype)

    %pad the length of the signal to be a multiple of M and a
    if mod(length(signal),lcm(a,M)) ~= 0
        signal = [signal; zeros(w-mod(length(signal),lcm(a,M)),1)];
    end

    % setup frame
    g = gabtight(wtype, a, M, w);
    Fr =  frame('dgtreal',g,a,M, phasetype);
    
    % get corrupted spectrogram
    coefs = frana(Fr, signal);
    spectrogram = framecoef2native(Fr,coefs); 
    
    s = 106; % gap start
    f = 111; % gap finish
    
    spectrogram(:,s:f) = 0;

    
    % calculate suboptimally restricted spectrogram
    N = size(spectrogram,2);
    [q,Q,p,P,S,F,u,v,U,V,L] = min_sgram_supp(w,a,M,s,f,N,phasetype);
    
    % setup fuctions for STFT and inverse STFT
    param.F = @(x) framecoef2native(Fr,frana(Fr, x));
    param.F_adj = @(z) frsyn(Fr,framenative2coef(Fr,z));


    subRestrict = spectrogram(:,q:Q);

    % settings for generalized Chambolle-Pock algorithm
    paramsolver.tau = 1;  % step size
    paramsolver.sigma = 1;  % step size
    paramsolver.eta = 1; % step size
    paramsolver.alpha = 1;  % relaxation parameter
    paramsolver.lambda = 0.1; % threshold (regularization parameter)
    paramsolver.I = 500; % number of iterations

    % inpaint the original and restricted sgram
    solution_orig = inpaint(spectrogram, param, paramsolver);
    solution_subR = inpaint(subRestrict, param, paramsolver);
    
    % calculate SNR in useful region
    ref = solution_orig(:,p:P);
    pred = solution_subR(:,U:V);
    
    SNR = snr(ref,ref-pred);
end