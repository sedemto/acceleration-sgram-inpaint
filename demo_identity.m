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

%% get shortened spectrogram
N = size(spectrogram,2);
[q,Q,p,P,S,F,u,v,U,V,L] = opt_sgram_supp(w,a,M,s,f,N,phasetype);

restricted = spectrogram(:,q:Q);

%% compare both sgrams after projection onto consistent sgram
sgram_sig = frsyn(Fr,framenative2coef(Fr,spectrogram));
projection_sgram = framecoef2native(Fr,frana(Fr,sgram_sig));

restriction_sig = frsyn(Fr,framenative2coef(Fr,restricted));
projection_restricted =  framecoef2native(Fr,frana(Fr,restriction_sig));

% calculate the SNR between the useful region of original and restricted sgram
ref = projection_sgram(:,p:P);
pred = projection_restricted(:,U:V);
SNR_useful_region = snr(ref, ref-pred);

disp("SNR in useful regions: "+SNR_useful_region);

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

% projection onto consistent sgram
new_restriction_sig = frsyn(Fr,framenative2coef(Fr,new_restriction));
projection_new_restricted =  framecoef2native(Fr,frana(Fr,new_restriction_sig));

% calculate the SNR using new useful region defined by shorter q
pred_new = projection_new_restricted(:,U_new:V_new);
SNR_shorter_restriction = snr(ref, ref-pred_new);

disp("SNR using shorter restriction: "+SNR_shorter_restriction);