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
[q,Q,p,P,S,F,u,v,U,V,L] = min_sgram_supp(w,a,M,s,f,N,phasetype);

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
solution_orig = inpaint(spectrogram, param, paramsolver);

% get reconstruction of restricted sgram
solution_restrict = inpaint(restricted, param, paramsolver);

%
ref = solution_orig(:,p:P);
pred = solution_restrict(:,U:V);
MSE = mean(abs((ref(:)-pred(:)).^2));

difference = ref-pred;
disp("MSE between the two reconstructions: "+MSE);

%% test shorter region
% test if q and Q is 1 shorter
q_shorter = q+1;
Q_shorter = P +w/a-2;

% move Q so that length of shortened sgram is multiple of M and a
if mod(Q_shorter-q_shorter+1,M/a) ~= 0
    Q_shorter = Q_shorter + M/a - mod(Q_shorter-q_shorter+1,M/a);
end

U_new = p-q_shorter+1;
V_new = U_new+P-p;


new_restriction = spectrogram(:,q_shorter:Q_shorter);

% get reconstruction of restricted sgram
solution_new_restrict = inpaint(new_restriction, param, paramsolver);

% calculate the MSE in the new useful region defined by shorter q
pred_new = solution_new_restrict(:,U_new:V_new);

MSE_shorter_restriction = mean(abs((ref(:)-pred_new(:)).^2));
diff2 = ref-pred_new;
disp("MSE of shorter restriction: "+MSE_shorter_restriction);