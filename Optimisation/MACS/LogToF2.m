function [y,ecc,a,The0] = LogToF2(x,r1,r2,chord,Dtheta,eccF,ic,eTmin,eTmax,pF)
    
    ExpX = exp((1/eTmin+1/eTmax)*x);
    eccT = eTmin*(ExpX-1)/(1+eTmin*ExpX/eTmax);

    ecc = sqrt(eccF^2+eccT^2);
    eccv = eccF*ic+eccT*[-ic(2); ic(1)];

    The0 = atan2(eccv(2),eccv(1));

    p = pF - r1*r2*eccT*sin(Dtheta)/chord;
    a = p/(1-ecc^2);
    ni1 = -The0;
    ni2 = Dtheta-The0;
    if ecc < 1-1.e-6
        E1 =EccAnomal(ni1,ecc);
        E2 =EccAnomal(ni2,ecc);
        while E2 < E1, E2 = E2+2*pi; end
        t12 = ToFell(E1,E2,a,ecc);
    elseif abs(ecc-1)<1.e-6
        E1 =EccAnomal(ni1,ecc-1.e-3);
        E2 =EccAnomal(ni2,ecc-1.e-3);
        while E2 < E1, E2 = E2+2*pi; end
        t12m = ToFell(E1,E2,a,ecc-1.e-3);
        F1 =HypAnomal(ni1,ecc+1.e-3);
        F2 =HypAnomal(ni2,ecc+1.e-3);
        while F2 < F1, F2 = F2+2*pi; end
        t12p = ToFhyp(F1,F2,a,ecc+1.e-3);
        t12 = (t12p+t12m)/2;        
    else
        F1 =HypAnomal(ni1,ecc);
        F2 =HypAnomal(ni2,ecc);
        while F2 < F1, F2 = F2+2*pi; end
        t12 = ToFhyp(F1,F2,a,ecc);
    end
    y = log(1*t12);
end
