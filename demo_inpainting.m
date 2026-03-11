clc;clear;close all

ltfatstart;
%% STFT parameters
w = 2048; % window length
a = w/4; % hop size
M = w; % number of freq. rows
wtype = 'hann'; % window type
phasetype = 'timeinv'; % freqinv or timeinv

% setup tight window
g = gabtight(wtype, a, M, w);

% STFT
Fr =  frame('dgtreal',g,a,M, phasetype);


%% create spectrogram and gap
[signal,fs] = audioread("a08_violin.wav");

%pad the length of the signal to be a multiple of M and a
if mod(length(signal),lcm(a,M)) ~= 0
    signal = [signal; zeros(w-mod(length(signal),lcm(a,M)),1)];
end

coefs = frana(Fr, signal);
spectrogram = framecoef2native(Fr,coefs); 

s = 106; % gap start
f = 111; % gap finish

spectrogram(:,s:f) = 0;

%% shortened spectrogram
N = size(spectrogram,2);
[q,Q,p,P,S,F,u,v,U,V,L] = opt_sgram_supp(w,a,M,s,f,N,phasetype);

restricted = spectrogram(:,q:Q);

%% spectrogram inpainting reconstruction comparison

% setup fuctions for STFT and inverse STFT
param.F = @(x) framecoef2native(Fr,frana(Fr, x));
param.F_adj = @(z) frsyn(Fr,framenative2coef(Fr,z));

% settings for generalized Chambolle-Pock algorithm
paramsolver.tau = 1;  % step size
paramsolver.sigma = 1;  % step size
paramsolver.eta = 1; % step size
paramsolver.alpha = 1;  % relaxation parameter
paramsolver.lambda = 0.1; % threshold (regularization parameter)
paramsolver.I = 500; % number of iterations


% get reconstruction of original sgram
tic
solution_orig = inpaint(spectrogram, param, paramsolver);
toc

% get reconstruction of restricted sgram
tic
solution_restrict = inpaint(restricted, param, paramsolver);
toc

%
ref = solution_orig(:,p:P);
pred = solution_restrict(:,U:V);
% MSE = mean(abs((ref(:)-pred(:)).^2));
SNR = snr(ref,pred-ref);

disp("SNR between the two reconstructions: "+SNR);
