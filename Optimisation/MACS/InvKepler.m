function  [ni2,err] = InvKepler(mu,a,ecc,ni1,T12)

err = 0;

nAVG = sqrt(mu/a^3);
ni2 = T12*3600*nAVG + ni1;

if ecc ~= 0
    Torb = (2*pi/nAVG)/3600;
    T2 = mod(T12,Torb);
    if ni1 ~= 0
        E1 = EccAnomal(ni1,ecc);
        T2 = T2 + (E1-ecc*sin(E1))/nAVG;
    end
    TorbMEZZI = Torb/2;
    E2k = EccAnomal(ni2,ecc);
    T2k = (E2k-ecc*sin(E2k))/nAVG/3600;
%    disp([E2k/2/pi T2k/Torb T2/Torb])
    k = 0;
    while abs(T2k-T2) > 1.e-12 && k < 1000        
        k = k+1;
        E2k = pi+(T2-TorbMEZZI)*(E2k-pi)/(T2k-TorbMEZZI);
        T2k = (E2k-ecc*sin(E2k))/nAVG/3600;
%    disp([E2k/2/pi T2k/Torb T2/Torb])
%    pause
    end
    if k >999
        disp('Non sono arrivato a convergenza in INVKEPLER')
        err = 1;
    end    
    ni2 = 2*atan2(sqrt((1+ecc)/(1-ecc))*sin(E2k/2),cos(E2k/2));
end

% disp([T12/2/pi,ni2/2/pi])
