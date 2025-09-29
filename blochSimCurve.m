function z_simulation = blochSimCurve(Freq, pools, pulseCell)
%   Change it to 36 for challenge 6&7
    nPulseRepeat = 1;
    % pulse: [strength in Hz, phase(0:x, pi/2:y, pi:-x), duration in s]
    magneticField = 3;       % in Tesla
    gyro = 1;                % gyromagnetic ratio relative to 1H



     Nfit=size(Freq,1);
     Nreadout=(size(pools,1)*2+1);
     watersignal = zeros(Nfit,2);
     
    for i=1:Nfit

        frequency=Freq(i);
        y1 = 0;


        y2Mat = blochSolve(pools, frequency, pulseCell, magneticField, nPulseRepeat, gyro, y1);
        % y2Mat: three dimensional matrix of magnetization vector after the pulse
        %     dim 1: [x0; x1; x2; ...y0; y1; y2;...z0; z1; z2; ...]
        %     dim 2: [after pulse 1; after pulse 2; after pulse 3;...]
        %     dim 3: [after repeat 1; after repeat 2; after repeat 3;...]
        % sequence: y1 -> [pulse1, pulse2,...] -> [pulse1, pulse2,...] -> ...
        %                 \_________________________________________________/
        %                              repeat pulseRepeat times
        %     The calculation is in a rotating frame with a frequency of the pulse
        %plot result
        watersignal(i,1)=frequency;
        %for three pools. x1 x2 x3 y1 y2 y3, z1, z2 z3. so index is 7 for water
        watersignal(i,2)=y2Mat( Nreadout,end,end);

    end
    
    z_simulation = watersignal(:,2);
end

