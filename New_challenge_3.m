%% Simulate z spectrum for Challenge 1

addpath('toolbox');
clear all;
clc;
close all;

% pool    = {name,            T1,      T2, lifetime, chemicalShift, concentration}
freewater  = {'freewater',              1.0,            0.040,              1,              0 ,                             1.11e5};
cest1        = {'amide',                    1.0,                0.1,         1/50,            3.5,     100e-3 / 111 * 1.11e5};
mtc          = {'MTC',                     1.0,          4.0e-5,          1/30,           -3.0,            15 / 111 * 1.11e5,   'Lorentzian'};
cest2        = {'guanidine_1',          1.0,                0.1,     1/1000,            2.0,     100e-3 / 111 * 1.11e5};
% cest3        = {'amine_1',                1.0,                0.1,     1/5000,            3.0,     100e-3 / 111 * 1.11e5};
% cest4        = {'OH_1',                   1.0,                0.1,     1/3500,            1.0,     100e-3 / 111 * 1.11e5};
noe          = {'NOE',                     1.3,            0.005,         1/20,           -3.0,     500e-3 / 111 * 1.11e5};

%always set the third one as MTC  
% pools = {freewater; cest1; mtc;  cest2; cest3; cest4; noe};
pools = {freewater; cest1; mtc;  cest2; noe};

%   Frequency offset (ppm)    
Freq_list = -15:0.1:15;
Freq_list = [-300, Freq_list];
Freq = Freq_list';

%   Pulse
pulse(:,1) = 2.0 *  267.5154109126/2/pi;  % power of the pulse is 2 uT
pulse(:,2) = 0;   % phase
pulse(:,3) = 2;   % duration
pulseCell = {pulse};


z_simulation = blochSimCurve(Freq, pools, pulseCell);
challenge_1_results(:, 1) = Freq;
challenge_1_results(:, 2) = z_simulation;
save challenge_3_results

figure();
plot(challenge_1_results(:, 1), challenge_1_results(:, 2));
xlim([-15, 15]);