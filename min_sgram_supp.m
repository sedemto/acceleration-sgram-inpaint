function [q,Q,p,P,S,F,u,v,U,V,L] = min_sgram_supp(w,a,M,s,f,N, phasetype)
% Find the minimum range of spectrogram that still carries all info about the gap
%
%   Usage: [q,Q] = min_sgram_supp(w,a,M,s,f,N,phasetype);
%          computes the first and last index of spectrogram needed for
%          reconstruction
%
%
%          [q,Q,p,P,S,F,u,v,U,V,L] = min_sgram_supp(w,a,M,s,f,N,phasetype);
%
%
%   Input parameters:
%       w     : Window length.
%       a     : Hop size.
%       M     : Number of frequency channels.
%       s     : Left edge index of the gap.
%       f     : Right edge index of the gap.
%       N     : Length of the original spectrogram.
%   phasetype : 'freqinv' or 'timeinv'
%
%         
%   Output parameters:
%       q     : First index of shortened spectrogram wrt. original spectrogram.
%       Q     : Last index of shortened spectrogram wrt. original spectrogram.
%       p     : Index of the first "useful" window overlapping with gap.
%       P     : Index of the last "useful" window overlapping with gap.
%       S     : Index of the first "useful" window within time-domain signal.
%       F     : Index of the last "useful" window within time-domain signal.
%       u     : Left edge of the gap within shortened spectrogram.
%       v     : Right edge of the gap within shortened spectrogram.
%       U     : Index of first "useful" window within shortened spectrogram.
%       V     : Index of last "useful" window within shortened spectrogram.
%       L     : Length of the shortened spectrogram.

% Author: Peter Balušík (2026)
% Brno University of Technology, Czech Republic

%%
% index of first useful window overlapping with the gap
p = s-w/a+1;

% index of last useful window overlapping with the gap
P = f+w/a-1;

% first index of the shortened part
q = p-w/a+1;

% last index of the shortened part
Q = P+w/a-1;

% if phasetype is 'freqinv' move q to correct position
if (phasetype=='freqinv')
    if mod(q-1,w/a) ~= 0
        q = q - mod(q-1,w/a);
    end
end


% move Q so that length of shortened sgram is multiple of M and a
if mod(Q-q+1,M/a) ~= 0
    Q = Q + M/a - mod(Q-q+1,M/a);
end


% checks if q<1 and Q>N
if q<1
    disp('Index q exceeds the length of the original spectrogram, q < 1!')
end

if Q>N
    disp('Index Q exceeds length of the original spectrogram!') %It can be padded by zeros
end

% position of the gap within the shortened spectrogram
u = s-q+1;
v = f-q+1;

% index of the first and last useful window within the shortened spectrogram
U = p-q+1;
V = U+P-p;

%final length of the spectrogram used
L = Q-q+1;

% calculate the first index of the shortened sgram within a time-domain signal
S = (q-1)*a +1;

% calculate the last index of the shortened sgram within a time-domain signal
F = Q*a;

end