function y2Mat = blochSolve(pools, frequency, pulseCell, magneticField, nPulseRepeat, gyro, y1)
% bloch simulation of the magnetization transfer 
% Author: Haifeng Zeng, hzeng6@jhmi.edu
% pool    = {name, T1, T2, lifetime, chemicalShift, concentration};
%    name: pool name
%    T1: T1 of the pool
%    T2: T2 of the pool
%    lifetime: lifetime of the exchangable proton in that pool, 
%              1/(exchange rate), 1 for free water 
%    chemicalShift: chemical shift of the exchangable proton, in ppm
%    concentration: concentration of the exchangable proton, in mM
% pools = {pool1; pool2; pool3; pool4;...}; pool1 is always free water pool
% frequency: carrier frequency of the pulse, in ppm
% pulseCell: {pulse1, pulse2, pulse3,...}
%     pulse: [strength1,phase1,duration1;strength2,phase2,duration2;...]
%         strength: in Hz 
%             0: delay, no pulse
%         phase: (0:x, pi/2:y, pi:-x, 3*pi/2:-y) 
%         duration: in s 
%             -1: crush transverse magnetization
%             inf: return to thermal equlibrium magnetization
% magneticField: magnetic field in Tesla
% nPulseRepeat: number of repeat of the pulses in pulseCell
% gyro: relative gyromagnetic ratio to 1H
% y1: initial magnetization vector in the follwoing order
%     [x(freewater); x(exchangable1); x(exchangable2);...
%      y(freewater); y(exchangable1); y(exchangable2);...
%      z(freewater); z(exchangable1); z(exchangable2);...]
%      z(freewater) = 1;
%     0: thermal equilibrium
% y2Mat: three dimensional matrix of magnetization vector after the pulse
%     dim 1: [x0; x1; x2; ...y0; y1; y2;...z0; z1; z2; ...]
%     dim 2: [after pulse 1; after pulse 2; after pulse 3;...]
%     dim 3: [after repeat 1; after repeat 2; after repeat 3;...]
% sequence: y1 -> [pulse1, pulse2,...] -> [pulse1, pulse2,...] -> ...
%                 \_________________________________________________/
%                              repeat pulseRepeat times
%     The calculation is in a rotating frame with a frequency of the pulse

gyro1H = 267.5154109126/2/pi;  % gyromagnetic ratio of 1H, MHz/T, 42.576400s 
pools = pools(:);
freewater = pools{1};
otherPools = pools(2:end);
nPools = length(pools);
nPulses = length(pulseCell);

if nPools > 2
        MTC_pool = pools{3, 1};
        mtc_lineshape = MTC_pool(7);    % the lineshape of the MTC background
end

%  Calculate the initial pool concentrations
%  concList: array of concentrations of all pools
concList = cellfun(@(pool) pool{6},  pools);
% y0: array of initial magnetizations, 
%     ([x0; x1; x2; ...y0; y1; y2;...z0; z1; z2; ...])
y0 = [zeros(nPools*2,1);concList(:)] / concList(1);
if (size(y1(:)) == size(y0))
    y1 = y1(:);
else
    y1 = y0;
end

% Calculate the exchange matrix
other2water = cellfun(@(pool) 1/pool{4},  otherPools);
water2other = cellfun(@(pool) pool{6}/pool{4},  otherPools) / freewater{6};
exchMat = zeros(nPools*3);
blockExchMat = zeros(nPools);
blockExchMat(1,2:end) = other2water(:)';
blockExchMat(2:end,1) = water2other(:);
for kPulCell = 0:2
    exchMat(kPulCell*nPools+1:(kPulCell+1)*nPools,...
        kPulCell*nPools+1:(kPulCell+1)*nPools) = blockExchMat;
end
exchMat = exchMat - diag(sum(exchMat,1));

% exchMat(3,:)=0;
% exchMat(:,3)=0;
% 
% exchMat(6,:)=0;
% exchMat(:,6)=0;

% Calculate the relaxation matrix
r1List = cellfun(@(pool) 1/pool{2},  pools);
r2List = cellfun(@(pool) 1/pool{3},  pools);
rMat = -diag([r2List;r2List;r1List]);

%
% rMat(3,3)=-1;
% rMat(6,6)=-1;
%

% Calculate the chemical shift matrix
% chemicalShiftList in ppm
chemicalShiftList = cellfun(@(pool) pool{5},  pools); 
% wList: chemical shift evolution frequency in radian
wList = 2*pi*magneticField*gyro1H*gyro*(chemicalShiftList-frequency);        
chemicalShiftMat = diag([-wList;zeros(nPools,1)],nPools) + ...
                   diag([wList;zeros(nPools,1)],-nPools);
            
%
% chemicalShiftMat(:,3)=0;
% chemicalShiftMat(3,:)=0;
% chemicalShiftMat(:,6)=0;
% chemicalShiftMat(6,:)=0;          
%

% Calculate the evolution matrix for the pulse
% http://en.wikipedia.org/wiki/Rotation_matrix
% section "9.3 Exponential map"
ypListCell = cell(1,nPulses);
transMatListCell = cell(1,nPulses);
crusherMat = diag([zeros(nPools*2,1);ones(nPools,1)]);
for kPulCell = 1:nPulses
    pulse = pulseCell{kPulCell};
    for kPul = 1:size(pulse,1)
        xList = ones(nPools,1) * pulse(kPul,1) * cos(pulse(kPul,2)) * pi * 2;
        yList = ones(nPools,1) * pulse(kPul,1) * sin(pulse(kPul,2)) * pi * 2;
        zList = zeros(nPools,1);
        pulseMat = diag([-zList;-xList], nPools) + diag(yList, nPools*2) + ...
                   diag([zList;xList], -nPools) + diag(-yList, -nPools*2);
               
        %
%         pulseMat(:,9)=0;
%         pulseMat(9,:)=0;
        %
            
         %first parameter R2mtc, second is the saturation power in Hz, third one is the
        %offset in Hz
        if pulse(kPul,1) == 0
            noRF_chemicalShiftList = cellfun(@(pool) pool{5},  pools); 
            % wList: chemical shift evolution frequency in radian
            wList = 2*pi*magneticField*gyro1H*gyro*(noRF_chemicalShiftList+1e6);        
            noRF_chemicalShiftList = diag([-wList;zeros(nPools,1)],nPools) + ...
                                        diag([wList;zeros(nPools,1)],-nPools);
            evolMat = exchMat + rMat + noRF_chemicalShiftList + pulseMat;
        else
            evolMat = exchMat + rMat + chemicalShiftMat + pulseMat;
        end

        if exist('mtc_lineshape' , 'var') == 1
                mtcpulse=sqrt(xList(3)^2+yList(3)^2);
                mtc_rf=RF_MT(1/r2List(3),mtcpulse,wList(3), mtc_lineshape); 
                evolMat(9,9)= evolMat(9,9)-mtc_rf;
        end
    
        t = pulse(kPul,3);
        yp = zeros(size(y0));
        transMat = crusherMat;
       
        if (~isinf(t) && t~=-1)
            yp = mldivide(evolMat, -rMat*y0);
            transMat = expm(evolMat*t);
        else
            disp('Signal is crushed or reach equilibrium');
        end
        ypListCell{kPulCell}(:,kPul) = yp;
        transMatListCell{kPulCell}(:,:,kPul) = transMat;
    end
end

% Calculate the evolution matrix
% dy/dt = (exchMat + chemicalShiftMat + pulseMat)y + rMat(y-y0)
% y is the vector of magnetization, y0 is the vector of magnetization at
%     thermal equilibrium
% change the equation into the following form
% d(y+yp) = (exchMat + chemicalShiftMat + pulseMat + rMat)*(y+yp)
%         = evolMat*(y+yp)
%    evolMat = exchMat + chemicalShiftMat + pulseMat + rMat;
%    yp      = mldivide(evolMat, -rMat*y0);
% the solution for y is 
% y = expm(evolMat*t)*(y+yp) - yp

% Calculae the evolution of the magnetization
y2Mat = zeros(nPools*3, nPulses, nPulseRepeat);
y2 = y1(:);
for kRepeat = 1:nPulseRepeat
    for kPulCell = 1:length(pulseCell) % 1
        pulse = pulseCell{kPulCell};
        for kPul = 1:size(pulse,1) % 80
            % if kPul == size(pulse,1) && kRepeat == nPulseRepeat
            %     break
            % else
            t = pulseCell{kPulCell}(kPul,3);
            if isinf(t)
                y2=y0;
                disp('t=inf');
            elseif (t==-1)
                y2 = crusherMat * y2;
                 disp('t=0');
            else
                yp = ypListCell{kPulCell}(:,kPul);
                transMat = transMatListCell{kPulCell}(:,:,kPul);
                y2 = transMat*(y2+yp) - yp;
            end
            % end
        end
        y2Mat(:,kPulCell,kRepeat) = y2;
    end
end

end