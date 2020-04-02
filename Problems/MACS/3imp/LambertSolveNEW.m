% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [vr1,vt1,vr2,vt2,ecc,a,p,The0,icont,ierr,t12err] = LambertSolveNEW(Dtheta,rho,t12des)

%disp([Dtheta,rho,t12des]);
    
ierr = 0;

r1 = 1;
r2 = rho;

rv1 = [r1;0];
rv2 = [r2*cos(Dtheta);r2*sin(Dtheta)];

% chord length
chord = norm(rv2-rv1);
% chord unit vector
ic = (rv2-rv1)/chord;

% Minimum eccentricity orbit (Min.Ecc.Orb.)
eccF = (r1-r2)/chord;
eccFv = eccF*ic;
pF = (r1+r2+dot(eccFv,rv1+rv2))/2;

a = pF/(1-eccF^2);
ecc = abs(eccF);

% Perifocal frame of Min.Ecc.Orb.
The0F = atan2(eccFv(2),eccFv(1));
ni1 = -The0F;
ni2 = Dtheta-The0F;
The0 = The0F;

% T.o.F. from 1 to 2 along Min.Ecc.Orb.
E1 =EccAnomal(ni1,ecc);
E2 =EccAnomal(ni2,ecc);
while E2 < E1, E2 = E2+2*pi; end
t12 = ToFell(E1,E2,a,ecc);
% Initial value of the transverse eccentricity vector component
x0 = 0;
y0 = log(t12);
ydes = log(1*t12des);

% Limits of the transverse eccentricity vector component
eTmax = sqrt(1-eccF^2);

% Newton-Raphson parameters

eps = 1.e-3;
tol = 1.e-10*10^(-log10(t12des)-log10(eTmax)-log10(min([Dtheta 2*pi-Dtheta])));
imax = 200;

% Solution for Dtheta > pi
if Dtheta > pi
    eccMAX = 1/cos(Dtheta/2);
    eTmin = sqrt(eccMAX^2-eccF^2);
    % Newton--Raphson Method
    icont = 0;
    [y0,ecc,a,p,The0] = LogToF2(x0,r1,r2,chord,Dtheta,eccF,ic,eTmin,eTmax,pF);
    while abs(y0-ydes) > tol && icont < imax
        icont = icont + 1;
        x1 = x0 + eps;
        if ~isreal(x1) 
           keyboard 
        end
        
        [y1,ecc,a,p,The0] = LogToF2(x1,r1,r2,chord,Dtheta,eccF,ic,eTmin,eTmax,pF);
        dydx = (y1-y0)/(x1-x0);
%       disp([x1,y1,dydx])
%         if dydx == 0
%             disp([Dtheta,rho,t12des])
%         end
        x0 = x0 + (ydes-y0)/dydx;
        if ~isreal(x0) 
           keyboard 
        end
        [y0,ecc,a,p,The0] = LogToF2(x0,r1,r2,chord,Dtheta,eccF,ic,eTmin,eTmax,pF);
%         disp([y0 ydes icont])
%         disp([abs(y0-ydes) > tol , icont < imax])
%         disp(abs(y0-ydes) > tol && icont < imax)
    end
%     disp(abs(y0-ydes) > tol && icont < imax)
%     disp('==========')
else
    % Newton--Raphson Method
    icont = 0;
    [y0,ecc,a,p,The0] = LogToF1(x0,r1,r2,chord,Dtheta,eccF,ic,eTmax,pF);
    while abs(y0-ydes) > tol && icont > -imax
        icont = icont - 1;
%         x1 = x0*(1 + eps);
        x1 = x0 + eps;
        [y1,ecc,a,p,The0] = LogToF1(x1,r1,r2,chord,Dtheta,eccF,ic,eTmax,pF);
        dydx = (y1-y0)/(x1-x0);
%         if dydx == 0
%             disp([Dtheta,rho,t12des])
%         end
        x0 = x0 + (ydes-y0)/dydx;
        [y0,ecc,a,p,The0] = LogToF1(x0,r1,r2,chord,Dtheta,eccF,ic,eTmax,pF);
    end
end

if abs(icont) > imax-1
    disp('Non sono arrivato a convergenza')
%     disp([Dtheta,rho,t12des])
    ierr = 1;
    vr1 = 0;
    vt1 = 0;
    vr2 = 0;
    vt2 = 0;
else
    vr1 = ecc*sin(-The0)/sqrt(p);
    vt1 = sqrt(p);
    vr2 = ecc*sin(Dtheta-The0)/sqrt(p);
    vt2 = vt1/r2;
end
t12err = abs(exp(y0)-t12des);

end

%%%%%%%%%%
% Case 1 % Dtheta <= pi
%%%%%%%%%%

function [y,ecc,a,p,The0] = LogToF1(x,r1,r2,chord,Dtheta,eccF,ic,eTmax,pF)
    
    eccT = eTmax*(1-exp(-x/eTmax));
    ecc = sqrt(eccF^2+eccT^2);
    eccv = eccF*ic+eccT*[-ic(2); ic(1)];
    
    The0 = atan2(eccv(2),eccv(1));

    p = max(pF - r1*r2*eccT*sin(Dtheta)/chord,0);
    ni1 = -The0;
    ni2 = Dtheta-The0;
    if ecc < 1%1-1.e-8
        a = p/(1-ecc^2);
        E1 =EccAnomal(ni1,ecc);
        E2 =EccAnomal(ni2,ecc);
        while E2 < E1, E2 = E2+2*pi; end
        t12 = ToFell(E1,E2,a,ecc);
    elseif ecc == 1%abs(ecc-1)<1.e-8
        E1 =EccAnomal(ni1,ecc);%E1 =EccAnomal(ni1,ecc-1.e-8);
        E2 =EccAnomal(ni2,ecc);%E2 =EccAnomal(ni2,ecc-1.e-8);
        while E2 < E1, E2 = E2+2*pi; end
        a = p/(1-(ecc)^2);%a = p/(1-(ecc-1.e-8)^2);
        t12m = ToFell(E1,E2,a,ecc);%t12m = ToFell(E1,E2,a,ecc-1.e-8);
        F1 =HypAnomal(ni1,ecc);%F1 =HypAnomal(ni1,ecc+1.e-8);
        F2 =HypAnomal(ni2,ecc);%F2 =HypAnomal(ni2,ecc+1.e-8);
        while F2 < F1, F2 = F2+2*pi; end
        a = p/(1-(ecc)^2);%a = p/(1-(ecc+1.e-8)^2);
        t12p = ToFhyp(F1,F2,a,ecc);%t12p = ToFhyp(F1,F2,a,ecc+1.e-8);
        t12 = (t12p+t12m)/2;
    else
        a = p/(1-ecc^2);
        F1 =HypAnomal(ni1,ecc);
        F2 =HypAnomal(ni2,ecc);
        while F2 < F1, F2 = F2+2*pi; end
        t12 = ToFhyp(F1,F2,a,ecc);
    end
    y = log(1*t12);
    
    if ~isreal(y)
       keyboard 
    end
    
%     disp([a,p,ecc])
%     disp([eccT,t12])
end

%%%%%%%%%%
% Case 2 % Dtheta > pi
%%%%%%%%%%

function [y,ecc,a,p,The0] = LogToF2(x,r1,r2,chord,Dtheta,eccF,ic,eTmin,eTmax,pF)
    
%     if x > 0
%         eccT = eTmax*(1-exp(-x/eTmax));
%     else
%         eccT = eTmin*(exp(x/eTmin)-1);
%     end

    ExpX = exp((1/eTmin+1/eTmax)*x);
    eccT = eTmin*(ExpX-1)/(1+eTmin*ExpX/eTmax);

    ecc = sqrt(eccF^2+eccT^2);
    eccv = eccF*ic+eccT*[-ic(2); ic(1)];

    The0 = atan2(eccv(2),eccv(1));

    p = max(pF - r1*r2*eccT*sin(Dtheta)/chord,0);
    ni1 = -The0;
    ni2 = Dtheta-The0;
    if ecc < 1%1-1.e-8
        a = p/(1-ecc^2);
        E1 =EccAnomal(ni1,ecc);
        E2 =EccAnomal(ni2,ecc);
        while E2 < E1, E2 = E2+2*pi; end
        t12 = ToFell(E1,E2,a,ecc);
    elseif ecc==1 %abs(ecc-1)<1.e-8
        E1 =EccAnomal(ni1,ecc);%E1 =EccAnomal(ni1,ecc-1.e-8);
        E2 =EccAnomal(ni2,ecc);%E2 =EccAnomal(ni2,ecc-1.e-8);
        while E2 < E1, E2 = E2+2*pi; end
        a = p/(1-(ecc)^2);%a = p/(1-(ecc-1.e-8)^2);
        t12m = ToFell(E1,E2,a,ecc);%t12m = ToFell(E1,E2,a,ecc-1.e-8);
        F1 =HypAnomal(ni1,ecc);%F1 =HypAnomal(ni1,ecc+1.e-8);
        F2 =HypAnomal(ni2,ecc);%F2 =HypAnomal(ni2,ecc+1.e-8);
        while F2 < F1, F2 = F2+2*pi; end
        a = p/(1-(ecc)^2);%a = p/(1-(ecc+1.e-8)^2);
        t12p = ToFhyp(F1,F2,a,ecc);%t12p = ToFhyp(F1,F2,a,ecc+1.e-8);
        t12 = (t12p+t12m)/2;        
    else
        a = p/(1-ecc^2);
        F1 =HypAnomal(ni1,ecc);
        F2 =HypAnomal(ni2,ecc);
        while F2 < F1, F2 = F2+2*pi; end
        t12 = ToFhyp(F1,F2,a,ecc);
    end
    y = log(1*t12);
    
    if ~isreal(y)
       keyboard 
    end

%     disp([a,p,ecc])
%     disp([eccT,t12])
end
