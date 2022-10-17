%% Simulate z spectrum for Challenge 2: two pools and 2 uT, 2 s

addpath('toolbox');
clear all;
clc;
close all;

% pool       = {name,                      T1,             T2,   lifetime, chemicalShift,                  concentration}
freewater   = {'freewater',              3.0,            2.0,              1,                   0,                             1.11e5};
cest1          = {'creatine_gua',      1.05,            0.1,         1/50,                 1.9,     55.5e-3 / 111 * 1.11e5};


%always set the third one as MTC  
pools = {freewater; cest1};

%   Frequency offset (ppm)    
Freq_list = -15:0.1:15;
Freq_list = [-300, Freq_list];
Freq = Freq_list';

%   Pulse
pulse(:,1) = 2.0 *  267.5154109126/2/pi;  % power of the pulse is 2 uT
pulse(:,2) = 0;   % phase
pulse(:,3) = 15;   % duration
pulseCell = {pulse};


z_simulation = blochSimCurve(Freq, pools, pulseCell);
challenge_2_results(:, 1) = Freq;
challenge_2_results(:, 2) = z_simulation;
save New_challenge_1_results

figure();
plot(challenge_2_results(:, 1), challenge_2_results(:, 2));
xlim([-15, 15]);