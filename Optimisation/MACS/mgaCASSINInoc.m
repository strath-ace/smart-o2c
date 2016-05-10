function [DVtot,rp] = mgaCASSINInoc(t)

%**************************************************************************
%Definition of the gravitational constants of the various planets
%(used in the powered swing-by routine) and of the sun (used in the lambert
%solver routine)
%**************************************************************************


mu(2)=324860;         %Gravitational constant of Venus
mu(3)=398601.19;      %Gravitational constant of Earth
mu(5)=126.7e6;        %Gravitational constant of Jupiter
mu(6)=37.9e6;         %Gravitational constant of Saturn
muSUN=1.32712428e+11; %Gravitational constant of Sun

%**************************************************************************
%Evaluation of position and velocities of the planets
%**************************************************************************

[r1,v1]=pleph_an ( real(t(1)) , 3 );
[r2,v2]=pleph_an ( real(t(1) + t(2)) , 2 );
[r3,v3]=pleph_an ( real(t(1) + t(2) + t(3)) , 2 );
[r4,v4]=pleph_an ( real(t(1) + t(2) + t(3) + t(4)) , 3 );
[r5,v5]=pleph_an ( real(t(1) + t(2) + t(3) + t(4) + t(5)) , 5 );
[r6,v6]=pleph_an ( real(t(1) + t(2) + t(3) + t(4) + t(5) + t(6)) , 6 );

%**************************************************************************
%Evaluation of the first Lambert arc and of the departure velocity
%**************************************************************************
 lw=vett(r1,r2);
    lw=sign(lw(3));
    if lw==1
        lw=0;
    else
        lw=1;
    end
    [V1,V21]=lambertI(r1,r2,t(2)*24*60*60,muSUN,lw);
    DVdep=norm(V1-v1);
    
%**************************************************************************
%Evaluation of the second Lambert arc and of the first Venus swingby DV
%**************************************************************************
 lw=vett(r2,r3);
    lw=sign(lw(3));
    if lw==1
        lw=0;
    else
        lw=1;
    end
    [V22,V31]=lambertI(r2,r3,t(3)*24*60*60,muSUN,lw); 
    Vin=norm(V21-v2);
    Vout=norm(V22-v2);
    alpha=acos((V21-v2)'*(V22-v2)/Vin/Vout);
    [DV1 , rp(1)]=PowSwingByInv(Vin,Vout,alpha);
    rp(1)=rp(1)*mu(2);
    
%**************************************************************************
%Evaluation of the third Lambert arc and of the second Venus swingby DV
%**************************************************************************
 lw=vett(r3,r4);
    lw=sign(lw(3));
    if lw==1
        lw=0;
    else
        lw=1;
    end
    [V32,V41]=lambertI(r3,r4,t(4)*24*60*60,muSUN,lw); 
    Vin=norm(V31-v3);
    Vout=norm(V32-v3);
    alpha=acos((V31-v3)'*(V32-v3)/Vin/Vout);
    [DV2 , rp(2)]=PowSwingByInv(Vin,Vout,alpha);
    rp(2)=rp(2)*mu(2);

%**************************************************************************
%Evaluation of the fourth Lambert arc and of the Earth swingby DV
%**************************************************************************
 lw=vett(r4,r5);
    lw=sign(lw(3));
    if lw==1
        lw=0;
    else
        lw=1;
    end
    [V42,V51]=lambertI(r4,r5,t(5)*24*60*60,muSUN,lw); 
    Vin=norm(V41-v4);
    Vout=norm(V42-v4);
    alpha=acos((V41-v4)'*(V42-v4)/Vin/Vout);
    [DV3 , rp(3)]=PowSwingByInv(Vin,Vout,alpha);
    rp(3)=rp(3)*mu(3);
    
%**************************************************************************
%Evaluation of the sixth Lambert arc and of the Jupiter swingby DV
%**************************************************************************
 lw=vett(r5,r6);
    lw=sign(lw(3));
    if lw==1
        lw=0;
    else
        lw=1;
    end
    [V52,V61]=lambertI(r5,r6,t(6)*24*60*60,muSUN,lw); 
    Vin=norm(V51-v5);
    Vout=norm(V52-v5);
    alpha=acos((V51-v5)'*(V52-v5)/Vin/Vout);
    [DV4 , rp(4)]=PowSwingByInv(Vin,Vout,alpha);
    rp(4)=rp(4)*mu(5);
    
%**************************************************************************
%Evaluation of the arrival DV 
%**************************************************************************

periSaturn=108950;  %Periapsis of the target orbit around Saturn 
                    %taken from the Cassini mission
e=.98;              %eccentricity of the target orbit

DVarr=norm(V61-v6); %Relative velocity at Saturn
DVper=sqrt(DVarr^2+2*mu(6)/periSaturn);  %Hyperbola
DVper2=sqrt(2*mu(6)/periSaturn-mu(6)/periSaturn*(1-e)); %Ellipse
DVarr=abs(DVper-DVper2);

DVtot=DVdep+DV1+DV2+DV3+DV4+DVarr;






%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************
%                       ASTRODYNAMICAL ROUTINES
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************

%Programmed by: The original version of this file was written in Fortran 77
%and is part of ESOC library. The file has later been translated under the
%name uplanet.m in Matlab, but that version was affected by several
%mistakes. In particular Pluto orbit was affected by great uncertainties.
%As a consequence its analytical eph have been substituted by a fifth order
%polynomial least square fit generated by Dario Izzo (ESA ACT). JPL405
%ephemerides (Charon-Pluto barycenter) have been used to produce the
%coefficients. WARNING: Pluto ephemerides should not be used outside the
%range 2000-2100;
%
%The analytical ephemerides of the nine planets of our solar system are
%returned in rectangular coordinates referred to the ecliptic J2000
%reference frame. Using analytical ephemerides instead of real ones (for
%examples JPL ephemerides) is faster but not as accurate, especially for
%the outer planets
%
%Date:                  28/05/2004  
%Revision:              1
%Tested by:             ----------
%
%
%Usage: [r,v,E]=pleph_an ( mjd2000 , IBODY );
%
%Inputs:
%           mjd2000:    Days elapsed from 1st of January 2000 counted from
%           midnight. 
%           planet:     Integer form 1 to 9 containing the number
%           of the planet (Mercury, Venus, Earth, Mars, Jupiter, Saturn,
%           Uranus, Neptune, Pluto) 
%
%Outputs:
%           r:         Column vector containing (in km) the position of the
%                      planet in the ecliptic J2000
%           v:         Column vector containing (in km/sec.) the velocity of
%                      the planet in the ecliptic J2000
%           E:         Column Vectors containing the six keplerian parameters,
%                      (a,e,i,OM,om,Eccentric Anomaly)
%
%Xref:      the routine needs a variable mu on the workspace whose 10th
%           element should be the sun gravitational parameter in KM^3/SEC^2


function  [r,v,E]=pleph_an ( mjd2000, planet);

RAD=pi/180;
AU  = 149597870.66;        %Astronomical Unit
KM=AU;


%
%  T = JULIAN CENTURIES SINCE 1900
%
T   = (mjd2000 + 36525.00)/36525.00;
TT  = T*T;
TTT = T*TT;
%
%  CLASSICAL PLANETARY ELEMENTS ESTIMATION IN MEAN ECLIPTIC OF DATE
%
switch planet
    %
    %  MERCURY
    %
    case 1 
        E(1) = 0.38709860;
        E(2) = 0.205614210 + 0.000020460*T - 0.000000030*TT;
        E(3) = 7.002880555555555560 + 1.86083333333333333e-3*T - 1.83333333333333333e-5*TT;
        E(4) = 4.71459444444444444e+1 + 1.185208333333333330*T + 1.73888888888888889e-4*TT;
        E(5) = 2.87537527777777778e+1 + 3.70280555555555556e-1*T +1.20833333333333333e-4*TT;
        XM   = 1.49472515288888889e+5 + 6.38888888888888889e-6*T;
        E(6) = 1.02279380555555556e2 + XM*T;
        %
        %  VENUS
        %
    case 2
        E(1) = 0.72333160;
        E(2) = 0.006820690 - 0.000047740*T + 0.0000000910*TT;
        E(3) = 3.393630555555555560 + 1.00583333333333333e-3*T - 9.72222222222222222e-7*TT;
        E(4) = 7.57796472222222222e+1 + 8.9985e-1*T + 4.1e-4*TT;
        E(5) = 5.43841861111111111e+1 + 5.08186111111111111e-1*T -1.38638888888888889e-3*TT;
        XM   = 5.8517803875e+4 + 1.28605555555555556e-3*T;
        E(6) = 2.12603219444444444e2 + XM*T;
        %
        %  EARTH
        %
    case 3
        E(1) = 1.000000230;
        E(2) = 0.016751040 - 0.000041800*T - 0.0000001260*TT;
        E(3) = 0.00;
        E(4) = 0.00;
        E(5) = 1.01220833333333333e+2 + 1.7191750*T + 4.52777777777777778e-4*TT + 3.33333333333333333e-6*TTT;
        XM   = 3.599904975e+4 - 1.50277777777777778e-4*T - 3.33333333333333333e-6*TT;
        E(6) = 3.58475844444444444e2 + XM*T;
        %
        %  MARS
        %
    case 4
        E(1) = 1.5236883990;
        E(2) = 0.093312900 + 0.0000920640*T - 0.0000000770*TT;
        E(3) = 1.850333333333333330 - 6.75e-4*T + 1.26111111111111111e-5*TT;
        E(4) = 4.87864416666666667e+1 + 7.70991666666666667e-1*T - 1.38888888888888889e-6*TT - 5.33333333333333333e-6*TTT;
        E(5) = 2.85431761111111111e+2 + 1.069766666666666670*T +  1.3125e-4*TT + 4.13888888888888889e-6*TTT;
        XM   = 1.91398585e+4 + 1.80805555555555556e-4*T + 1.19444444444444444e-6*TT;
        E(6) = 3.19529425e2 + XM*T;
        
        %
        %  JUPITER
        %
    case 5
        E(1) = 5.2025610;
        E(2) = 0.048334750 + 0.000164180*T  - 0.00000046760*TT -0.00000000170*TTT;
        E(3) = 1.308736111111111110 - 5.69611111111111111e-3*T +  3.88888888888888889e-6*TT;
        E(4) = 9.94433861111111111e+1 + 1.010530*T + 3.52222222222222222e-4*TT - 8.51111111111111111e-6*TTT;
        E(5) = 2.73277541666666667e+2 + 5.99431666666666667e-1*T + 7.0405e-4*TT + 5.07777777777777778e-6*TTT;
        XM   = 3.03469202388888889e+3 - 7.21588888888888889e-4*T + 1.78444444444444444e-6*TT;
        E(6) = 2.25328327777777778e2 + XM*T;
        %
        %  SATURN
        %
    case 6
        E(1) = 9.5547470;
        E(2) = 0.055892320 - 0.00034550*T - 0.0000007280*TT + 0.000000000740*TTT;
        E(3) = 2.492519444444444440 - 3.91888888888888889e-3*T - 1.54888888888888889e-5*TT + 4.44444444444444444e-8*TTT;
        E(4) = 1.12790388888888889e+2 + 8.73195138888888889e-1*T -1.52180555555555556e-4*TT - 5.30555555555555556e-6*TTT;
        E(5) = 3.38307772222222222e+2 + 1.085220694444444440*T + 9.78541666666666667e-4*TT + 9.91666666666666667e-6*TTT;
        XM   = 1.22155146777777778e+3 - 5.01819444444444444e-4*T - 5.19444444444444444e-6*TT;
        E(6) = 1.75466216666666667e2 + XM*T;
        %
        %  URANUS
        %
    case 7
        E(1) = 19.218140;
        E(2) = 0.04634440 - 0.000026580*T + 0.0000000770*TT;
        E(3) = 7.72463888888888889e-1 + 6.25277777777777778e-4*T + 3.95e-5*TT;
        E(4) = 7.34770972222222222e+1 + 4.98667777777777778e-1*T + 1.31166666666666667e-3*TT;
        E(5) = 9.80715527777777778e+1 + 9.85765e-1*T - 1.07447222222222222e-3*TT - 6.05555555555555556e-7*TTT;
        XM   = 4.28379113055555556e+2 + 7.88444444444444444e-5*T + 1.11111111111111111e-9*TT;
        E(6) = 7.26488194444444444e1 + XM*T;
        %
        %  NEPTUNE
        %
    case 8
        E(1) = 30.109570;
        E(2) = 0.008997040 + 0.0000063300*T - 0.0000000020*TT;
        E(3) = 1.779241666666666670 - 9.54361111111111111e-3*T - 9.11111111111111111e-6*TT;
        E(4) = 1.30681358333333333e+2 + 1.0989350*T + 2.49866666666666667e-4*TT - 4.71777777777777778e-6*TTT;
        E(5) = 2.76045966666666667e+2 + 3.25639444444444444e-1*T + 1.4095e-4*TT + 4.11333333333333333e-6*TTT;
        XM   = 2.18461339722222222e+2 - 7.03333333333333333e-5*T;
        E(6) = 3.77306694444444444e1 + XM*T;
        %
        %  PLUTO
        %
    case 9
        %Fifth order polynomial least square fit generated by Dario Izzo
        %(ESA ACT). JPL405 ephemerides (Charon-Pluto barycenter) have been used to produce the coefficients.
        %This approximation should not be used outside the range 2000-2100;
        T=mjd2000/36525;
        TT=T*T;
        TTT=TT*T;
        TTTT=TTT*T;
        TTTTT=TTTT*T;
        E(1)=39.34041961252520 + 4.33305138120726*T - 22.93749932403733*TT + 48.76336720791873*TTT - 45.52494862462379*TTTT + 15.55134951783384*TTTTT;
        E(2)=0.24617365396517 + 0.09198001742190*T - 0.57262288991447*TT + 1.39163022881098*TTT - 1.46948451587683*TTTT + 0.56164158721620*TTTTT;
        E(3)=17.16690003784702 - 0.49770248790479*T + 2.73751901890829*TT - 6.26973695197547*TTT + 6.36276927397430*TTTT - 2.37006911673031*TTTTT;
        E(4)=110.222019291707 + 1.551579150048*T - 9.701771291171*TT + 25.730756810615*TTT - 30.140401383522*TTTT + 12.796598193159 * TTTTT;
        E(5)=113.368933916592 + 9.436835192183*T - 35.762300003726*TT + 48.966118351549*TTT - 19.384576636609*TTTT - 3.362714022614 * TTTTT;
        E(6)=15.17008631634665 + 137.023166578486*T + 28.362805871736*TT - 29.677368415909*TTT - 3.585159909117*TTTT + 13.406844652829 * TTTTT;
end
%
%  CONVERSION OF AU INTO KM, DEG INTO RAD
%
E(1)     =     E(1)*KM;
for  I = 3: 6
    E(I)     = E(I)*RAD;
end
E(6)     = mod(E(6), 2*pi);

%Conversion from mean anomaly to eccentric anomaly via Kepler's equation
EccAnom=M2E(E(6),E(2));
E(6)=EccAnom;

%Calcolo velocit� e posizione nel sistema J2000
[r,v]=conversion(E);




function [r,v] = conversion (E)

% Parametri orbitali

muSUN=  1.327124280000000e+011;       %gravitational parameter for the sun

a=E(1);
e=E(2);
i=E(3);
omg=E(4);
omp=E(5);
EA=E(6);


% Grandezze definite nel piano dell'orbita

b=a*sqrt(1-e^2);
n=sqrt(muSUN/a^3);

xper=a*(cos(EA)-e);
yper=b*sin(EA);

xdotper=-(a*n*sin(EA))/(1-e*cos(EA));
ydotper=(b*n*cos(EA))/(1-e*cos(EA));

% Matrice di trasformazione da perifocale a ECI

R(1,1)=cos(omg)*cos(omp)-sin(omg)*sin(omp)*cos(i);
R(1,2)=-cos(omg)*sin(omp)-sin(omg)*cos(omp)*cos(i);
R(1,3)=sin(omg)*sin(i);
R(2,1)=sin(omg)*cos(omp)+cos(omg)*sin(omp)*cos(i);
R(2,2)=-sin(omg)*sin(omp)+cos(omg)*cos(omp)*cos(i);
R(2,3)=-cos(omg)*sin(i);
R(3,1)=sin(omp)*sin(i);
R(3,2)=cos(omp)*sin(i);
R(3,3)=cos(i);

% Posizione nel sistema inerziale

r=R*[xper;yper;0];
v=R*[xdotper;ydotper;0];


function E=M2E(M,e)
tol=1e-13;
err=1;
E=M+e*cos(M);   %initial guess
while err>tol
    Enew=E-(E-e*sin(E)-M)/(1-e*cos(E));
    err=abs(E-Enew);
    E=Enew;
end







%Programmed by:         Dario Izzo (Advanced Concept Team) Date:
%28/05/2004 Revision:              1 Tested by:             ----------
%
%
%This routine implements a new algorithm that solves Lambert's problem. The
%algorithm has two major characteristics that makes it favorable to other
%existing ones. 
%
%   1) It describes the generic orbit solution of the boundary condition
%   problem through the variable X=log(1+cos(alpha/2)). By doing so the
%   graphs of the time of flight become defined in the entire real axis and
%   resembles a straight line. Convergence is granted within few iterations
%   for all the possible geometries (except, of course, when the transfer
%   angle is zero). When multiple revolutions are considered the variable is
%   X=tan(cos(alpha/2)*pi/2).
%
%   2) Once the orbit has been determined in the plane, this routine
%   evaluates the velocity vectors at the two points in a way that is not
%   singular for the transfer angle approaching to pi (Lagrange coefficient
%   based methods are numerically not well suited for this purpose).
%
%   As a result Lambert's problem is solved (with multiple revolutions
%   being accounted for) with the same computational effort for all
%   possible geometries. The case of near 180 transfers is also solved
%   efficiently.
%
%   We note here that even when the transfer angle is exactly equal to pi
%   the algorithm does solve the problem in the plane (it finds X), but it
%   is not able to evaluate the plane in which the orbit lies. A solution 
%   to this would be to provide the direction of the plane containing the
%   transfer orbit from outside. This has not been implemented in this 
%   routine since such a direction would depend on which application the 
%   transfer is going to be used in.
%
%Usage: [v1,v2,a,p,theta,iter]=lambertI(r1,r2,t,mu,lw,N,branch)
%
%Inputs:
%           r1=Position vector at departure (column) 
%           r2=Position vector at arrival (column, same units as r1)
%           t=Transfer time (scalar)
%           mu=gravitational parameter (scalar, units have to be
%           consistent with r1,t units) 
%           lw=1 if long way is chosen    
%           branch='l' if the left branch is selected in a problem where N
%           is not 0 (multirevolution)
%           N=number of revolutions
%
%Outputs:
%           v1=Velocity at departure        (consistent units)
%           v2=Velocity at arrival
%           a=semi major axis of the solution
%           p=semi latus rectum of the solution 
%           theta=transfer angle in rad
%           iter=number of iteration made by the newton solver (usually 6)
%
%please report bugs to dario.izzo@esa.int



function [v1,v2,a,p,theta,iter]=lambertI(r1,r2,t,mu,lw,N,branch)

%Preliminary control on the function call
if nargin==5
    N=0;
end
if t<=0
    warning('Negative time as input')
    v1=NaN;
    v2=NaN;
    return
end


tol=1e-11;  %Increasing the tolerance does not bring any advantage as the 
%precision is usually greater anyway (due to the rectification of the tof
%graph) except near particular cases such as parabolas in which cases a
%lower precision allow for usual convergence.


%Non dimensional units
R=sqrt(r1'*r1);
V=sqrt(mu/R);
T=R/V;

%working with non-dimensional radii and time-of-flight
r1=r1/R;
r2=r2/R;
t=t/T;                     

%Evaluation of the relevant geometry parameters in non dimensional units
r2mod=sqrt(r2'*r2);
theta=real(acos((r1'*r2)/r2mod)); %the real command is useful when theta is very 
                                  %close to pi and the acos function could return complex numbers
if lw
    theta=2*pi-theta;
end
c=sqrt(1+r2mod^2-2*r2mod*cos(theta)); %non dimensional chord
s=(1+r2mod+c)/2;                      %non dimensional semi-perimeter
am=s/2;                               %minimum energy ellipse semi major axis
lambda=sqrt(r2mod)*cos(theta/2)/s;    %lambda parameter defined in BATTIN's book



%We start finding the log(x+1) value of the solution conic:
%%NO MULTI REV --> (1 SOL)
if N==0
    inn1=-.5233;    %first guess point
    inn2=.5233;     %second guess point
    x1=log(1+inn1);
    x2=log(1+inn2);
    y1=log(x2tof(inn1,s,c,lw,N))-log(t);
    y2=log(x2tof(inn2,s,c,lw,N))-log(t);
    
    %Newton iterations
    err=1;
    i=0;
    while ((err>tol) & (y1~=y2))
        i=i+1;
        xnew=(x1*y2-y1*x2)/(y2-y1);
        ynew=log(x2tof(exp(xnew)-1,s,c,lw,N))-log(t);
        x1=x2;
        y1=y2;
        x2=xnew;
        y2=ynew;
        err=abs(x1-xnew);
    end
    iter=i;
    x=exp(xnew)-1;
    
    
    %%MULTI REV --> (2 SOL) SEPARATING RIGHT AND LEFT BRANCH
else 
    if branch=='l'
        inn1=-.5234;
        inn2=-.2234;
    else
        inn1=.7234;
        inn2=.5234;
    end
    x1=tan(inn1*pi/2);
    x2=tan(inn2*pi/2);
    y1=x2tof(inn1,s,c,lw,N)-t;
    
    y2=x2tof(inn2,s,c,lw,N)-t;
    err=1;
    i=0;
    
    %Newton Iteration
    while ((err>tol) & (i<60) & (y1~=y2))
        i=i+1;
        xnew=(x1*y2-y1*x2)/(y2-y1);
        ynew=x2tof(atan(xnew)*2/pi,s,c,lw,N)-t;
        x1=x2;
        y1=y2;
        x2=xnew;
        y2=ynew;
        err=abs(x1-xnew);	    
    end
    x=atan(xnew)*2/pi;
    iter=i;
end

%The solution has been evaluated in terms of log(x+1) or tan(x*pi/2), we
%now need the conic. As for transfer angles near to pi the lagrange
%coefficient technique goes singular (dg approaches a zero/zero that is
%numerically bad) we here use a different technique for those cases. When
%the transfer angle is exactly equal to pi, then the ih unit vector is not
%determined. The remaining equations, though, are still valid.


a=am/(1-x^2);                       %solution semimajor axis
%calcolo psi
if x<1 %ellisse
    beta=2*asin(sqrt((s-c)/2/a));
    if lw
        beta=-beta;
    end
    alfa=2*acos(x);
    psi=(alfa-beta)/2;
    eta2=2*a*sin(psi)^2/s;
    eta=sqrt(eta2);
else %iperbole
    beta=2*asinh(sqrt((c-s)/2/a));
    if lw
        beta=-beta;
    end
    alfa=2*acosh(x);
    psi=(alfa-beta)/2;
    eta2=-2*a*sinh(psi)^2/s;
    eta=sqrt(eta2);
end
p=r2mod/am/eta2*sin(theta/2)^2;     %parameter of the solution
sigma1=1/eta/sqrt(am)*(2*lambda*am-(lambda+x*eta));
ih=vers(vett(r1,r2)');
if lw
    ih=-ih;
end

vr1 = sigma1;
vt1 = sqrt(p);
v1  = vr1 * r1   +   vt1 * vett(ih,r1)';

vt2=vt1/r2mod;
vr2=-vr1+(vt1-vt2)/tan(theta/2);
v2=vr2*r2/r2mod+vt2*vett(ih,r2/r2mod)';
v1=v1*V;
v2=v2*V;
a=a*R;
p=p*R;

%Subfunction that evaluates the time of flight as a function of x
function t=x2tof(x,s,c,lw,N)  

am=s/2;
a=am/(1-x^2);
if x<1 %ELLISSE
    beta=2*asin(sqrt((s-c)/2/a));
    if lw
        beta=-beta;
    end
    alfa=2*acos(x);
else   %IPERBOLE
    alfa=2*acosh(x);
    beta=2*asinh(sqrt((s-c)/(-2*a)));
    if lw
        beta=-beta;
    end
end
t=tofabn(a,alfa,beta,N);

%subfunction that evaluates the time of flight via Lagrange expression
function t=tofabn(sigma,alfa,beta,N)

if sigma>0
    t=sigma*sqrt(sigma)*((alfa-sin(alfa))-(beta-sin(beta))+N*2*pi);
else
    t=-sigma*sqrt(-sigma)*((sinh(alfa)-alfa)-(sinh(beta)-beta));
end

%subfunction that evaluates unit vectors
function v=vers(V) 
v=V/sqrt(V'*V);

%subfunction that evaluates vector product
function ansd = vett(r1,r2)  
ansd(1)=(r1(2)*r2(3)-r1(3)*r2(2));
ansd(2)=(r1(3)*r2(1)-r1(1)*r2(3));
ansd(3)=(r1(1)*r2(2)-r1(2)*r2(1));
  








%Programmed by: Dario Izzo (ESA/ACT)
%
%Date:                  24/09/2005  
%Revision:              1
%Tested by:             ----------
%
%
%Usage: [DV,rp,iter] = PowSwingByInv(Vin,Vout,alpha)
%
%Outputs:
%           DV:    Velcity Increment of the Powered SwingBy (non dimensional)
%           rp:    Pericenter radius found.
%           iter:  Number of iteration to converge (-1 if convergence is failed)
%
%Inputs:
%           Vin:   Incoming hyperbolic velocity modulus  (non dimensional)
%           Vout:  Outgoing hyperbolic velocity modulus  (non dimensional)
%           alpha: Angle between Vin and Vout (in rad.)
%
%Comments:  The non dimensional units are R for the length and sqrt(mu/R)
%for the velocity --> gravitational constant is one. R may be choosen
%freely as any relevant length. Magic of the non dimensional forms: if we
%forget about dimension and we call the routine, then DV is returned in the
%same units as the input parameters, and rp has to be multiplied by the
%planet gravitational constant (unit consistent with the velocity input)
%to be transformed in the length.


function [DV,rp,iter]=PowSwingByInv(Vin,Vout,alpha)


aIN=1/Vin^2;    %semimajor axis of the incoming hyperbola
aOUT=1/Vout^2;  %semimajor axis of the outcoming hyperbola

%We find the perigee radius with an iteration method based on the gradient
%of the function. Attention has to be payed to the initial point as the
%function is not defined for rp<0. The option here implemented considers
%halfing the perigee radius whenever the gradient pushes the next iteration
%in a non defined zone.
i=0;
maxiter=30;     %maximum number of iteration allowed for the gradient method
rp=1;           %Initial point
err=1;
stop=1;

while stop
    i=i+1;
    f=asin(aIN/(aIN+rp))+asin(aOUT/(aOUT+rp))-alpha;
    df=-aIN/sqrt(rp^2+2*aIN*rp)/(aIN+rp)-aOUT/sqrt(rp^2+2*aOUT*rp)/(aOUT+rp);
    rpNew=rp-f/df;
    rp=rpNew;
    err=abs(f);
    if rp<0
        rp=1e-8;
    end
    if err<1e-7|i>=maxiter|abs(df)<1e-8
        stop=0;
    end
end

%Evaluation of the DV
DV=abs(sqrt(Vout^2+2/rp)-sqrt(Vin^2+2/rp));

%If the maximum number of iteration is achieved the returned number of
%iteration is -1.
iter=i;
if iter==maxiter
    iter=-1;
end

