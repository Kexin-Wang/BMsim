%% Simulate z spectrum for Challenge 1

addpath('toolbox');
addpath('yaml_toolbox');
clear all;
clc;
close all;

yaml_file = '/Users/wangk6/Documents/postdoc/BMsimulation/case_4_5pool_model.yaml';
YamlStruct = ReadYaml(yaml_file);
%%
% pool    = {name,            T1,      T2, lifetime, chemicalShift, concentration}
freewater  = {'freewater',              YamlStruct.water_pool.t1,    YamlStruct.water_pool.t2,              1,              0 ,    YamlStruct.water_pool.f*1.11e5};
cest1        = {'amide',      YamlStruct.cest_pool.amide.t1,    YamlStruct.cest_pool.amide.t2,     1/YamlStruct.cest_pool.amide.k,       YamlStruct.cest_pool.amide.dw,     YamlStruct.cest_pool.amide.f * 1.11e5};
mtc          = {'MTC',       YamlStruct.cest_pool.MT.t1,    YamlStruct.cest_pool.MT.t2,         1/YamlStruct.cest_pool.MT.k,       YamlStruct.cest_pool.MT.dw,      YamlStruct.cest_pool.MT.f * 1.11e5,  'Lorentzian'};
cest2        = {'guanidine_1',  YamlStruct.cest_pool.guanidine_1.t1,    YamlStruct.cest_pool.guanidine_1.t2,     1/YamlStruct.cest_pool.guanidine_1.k,  YamlStruct.cest_pool.guanidine_1.dw,    YamlStruct.cest_pool.guanidine_1.f* 1.11e5};
% cest3        = {'amine_1',                1.0,                0.1,     1/5000,            3.0,     100e-3 / 111 * 1.11e5};
% cest4        = {'OH_1',                   1.0,                0.1,     1/3500,            1.0,     100e-3 / 111 * 1.11e5};
noe          = {'NOE',  YamlStruct.cest_pool.NOE_1.t1,    YamlStruct.cest_pool.NOE_1.t2,         1/YamlStruct.cest_pool.NOE_1.k,       YamlStruct.cest_pool.NOE_1.dw,     YamlStruct.cest_pool.NOE_1.f * 1.11e5};

%always set the third one as MTC  
% pools = {freewater; cest1; mtc;  cest2; cest3; cest4; noe};
pools = {freewater; cest1; mtc;  cest2; noe};

%   Frequency offset (ppm)    
Freq_list = -2:0.05:2;
Freq_list = [-300, Freq_list];
Freq = Freq_list';

%   Pulse
pulse(:,1) = 3.7 *  267.5154109126009/2/pi;  % power of the pulse is 2 uT
pulse(:,2) = 0;   % phase
pulse(:,3) = 5e-3;   % duration
pulse2(:,1) = 0;  % power of the pulse is 2 uT
pulse2(:,2) = 0;   % phase
pulse2(:,3) = -1;   % duration
pulse3(:,1) = 0;  % power of the pulse is 2 uT
pulse3(:,2) = 0;   % phase
pulse3(:,3) = 6.5e-3;   % duration
pulseCell = {pulse,pulse2,pulse3};


z_simulation = blochSimCurve(Freq, pools, pulseCell);
challenge_1_results(:, 1) = Freq;
challenge_1_results(:, 2) = z_simulation;
save challenge_4_results

figure();
plot(challenge_1_results(:, 1), challenge_1_results(:, 2));
xlim([-2, 2]);