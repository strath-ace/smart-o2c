function [DvTOT,eccTransf,rTransfMIN,datiplot,datirv] = DvTransf_3a(x,p)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2018 University of Strathclyde and Authors-----------

% ------------------------
%      Initialization
% ------------------------
invkc=0;
mu      = p(1);

% LEO orbit (start - Orbit #1)
rLEO    = p(2);
the0_C  = p(3);
wLEO    = p(4);

% High ecc. orbit (arrival - Orbit #2)
a2      = p(5);   % [km]
e2      = p(6);
p2      = p(7);       % [km]
the0_T  = p(8);
THEO    = p(9);

% Optimization parameters
Twait=x(1);   % waiting time on LEO orbit
Ttran1=x(2);   % transfer time to h.e. orbit
Dtheta1=x(3);   % transfer time to h.e. orbit
Rp=x(4);   % transfer time to h.e. orbit
Ttran2=x(5);   % transfer time to h.e. orbit
        
Ttot = Twait+Ttran1+Ttran2;
         
% Position of chaser at start point of transfer
r1 = rLEO;
the1 = mod(the0_C+wLEO*Twait*3600,2*pi);
         
% Arrival point - final
Tres_T = mod(Ttot,THEO);
[ni2_T,err] = InvKepler(mu,a2,e2,the0_T,Tres_T);
if err > 0
    disp('Non sono arrivato a convergenza in INVKEPLER')
    invkc=1;
end
r2 = p2/(1+e2*cos(ni2_T));
the2 = ni2_T; % Assumption: Orbit #2 periapsis reference for anomalies
datiplot=[r1 the1 Rp the1+Dtheta1 r2 the2];
% ------------------------
% Determine transfer orbit
% ------------------------

Dtheta2 = the2 - (the1+Dtheta1);      % transfer angle
while Dtheta2 < 0
    Dtheta2 = Dtheta2+2*pi;
end


epsilon=1.e-5;
if Dtheta2 < 2*pi-epsilon && Dtheta2 > epsilon

    % disp('--- 1 ---')
    % disp([Dtheta,r2/r1,Ttran*3600*wLEO])

    % Solve Lambert's problem
    [vr1a,vt1a,vr2a,vt2a,eccTransfa,aTransfa,pTransfa,ome0Transfa,iconta,ierra,t12errb] ...
        = LambertSolveNEW(Dtheta1,Rp/r1,Ttran1*3600*wLEO);
    if ~isreal([vr1a vt1a vr2a vt2a])
%         disp('Trovati valori di velocita'' con componente immaginaria')
%         disp([vr1a vt1a vr2a vt2a]')
        invkc=1;
    end
    if ierra > 0
%         disp('Non sono arrivato a convergenza in LAMBERT')
%         disp(x)
        invkc=1;
    end

    ome0Transfa = mod(ome0Transfa+the1,2*pi);
    if ome0Transfa>pi
        ome0Transfa=-2*pi+ome0Transfa;
    end
    % Transfer orbit minimum distance
    the1P = mod(the1,2*pi);
    the2P = mod(the1+Dtheta1,2*pi);

    rTransfMINa = rLEO;

    if the1P > the2P
        the1P = the1P - 2*pi;
    end
    if (the2P-ome0Transfa)*(ome0Transfa-the1P)>0
        rTransfMINa = aTransfa*rLEO*(1-eccTransfa);
    end

    Ttransf=2*pi*sqrt((Rp)^3/mu);
    wTransf=2*pi/Ttransf;
    
    vr1a = wLEO*rLEO*vr1a;
    vt1a = wLEO*rLEO*vt1a;
    vr2a = wLEO*rLEO*vr2a;
    vt2a = wLEO*rLEO*vt2a;

    vrI = 0;
    vtI = wLEO*rLEO;
    vrF = sqrt(mu/p2)*e2*sin(ni2_T);
    vtF = sqrt(mu/p2)*(1+e2*cos(ni2_T));
    if isnan(eccTransfa)
        ierra=2;
    end
    
    %**********************************************************************
    if ierra<=0
        [vr1b,vt1b,vr2b,vt2b,eccTransfb,aTransfb,pTransfb,ome0Transfb,icontb,ierrb,t12errb] ...
            = LambertSolveNEW(Dtheta2,r2/Rp,Ttran2*3600*wTransf);
        if ~isreal([vr1b vt1b vr2b vt2b])
%             disp('Trovati valori di velocita'' con componente immaginaria')
%             disp([vr1b vt1b vr2b vt2b]')
            invkc=1;
        end
        if ierrb > 0
%             disp('Non sono arrivato a convergenza in LAMBERT')
%             disp(x)
            invkc=1;
        end
        if isnan(eccTransfb)
            ierrb=2;
        end
        if ierrb<=0
            ome0Transfb = mod(ome0Transfb+the1+Dtheta1,2*pi);
            if ome0Transfb>pi
                ome0Transfb=-2*pi+ome0Transfb;
            end
            % Transfer orbit minimum distance
            the1P = mod(the1+Dtheta1,2*pi);
            the2P = mod(the2,2*pi);

            rTransfMINb = Rp;

            if the1P > the2P
                the1P = the1P - 2*pi;
            end
            if (the2P-ome0Transfb)*(ome0Transfb-the1P)>0
                rTransfMINb = aTransfb*Rp*(1-eccTransfb);
            end


            vr1b = wTransf*Rp*vr1b;
            vt1b = wTransf*Rp*vt1b;
            vr2b = wTransf*Rp*vr2b;
            vt2b = wTransf*Rp*vt2b;



            %DvTOTp = sqrt((vr1a-vrI)^2+(vt1a-vtI)^2)+sqrt((vrF-vr2b)^2+(vtF-vt2b)^2)+sqrt((vr1b-vr2a)^2+(vt1b-vt2a)^2);
            [DvTOT]=DVtot_3(the1,Dtheta1,the2,aTransfa*r1*1e3,eccTransfa,ome0Transfa,aTransfb*Rp*1e3,eccTransfb,ome0Transfb);
            DvTOT=DvTOT/1e3;
            datiplot=[datiplot eccTransfa,aTransfa,pTransfa,ome0Transfa  eccTransfb,aTransfb,pTransfb,ome0Transfb];
            datirv=[vr1a,vt1a,vr2a,vt2a vr1b,vt1b,vr2b,vt2b];
            rTransfMIN=min([rTransfMINa rTransfMINb]);
            eccTransf=max([eccTransfa eccTransfb]);
        else
            DvTOT = 1.e14;
            eccTransf = 10;
            rTransfMIN = 0;
        end
    else
        DvTOT = 1.e14;
        eccTransf = 10;
        rTransfMIN = 0;
    end
else
    DvTOT = 1.e14;
    eccTransf = 10;
    rTransfMIN = 0;
end
if invkc>0.5
    DvTOT = 1.e14;
    eccTransf = 10;
    rTransfMIN = 0;
end