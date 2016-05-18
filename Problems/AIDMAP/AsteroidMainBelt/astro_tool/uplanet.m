function  [E,ksun] = uplanet(mjd2000, IBODY)

%  ANALYTICAL EPHEMERIDES FOR PLANETS - P. DYSLI, 1977.
%  PLANETARY ORBITAL ELEMENTS ARE RESTITUTED IN A SUN-CENTRED (EARTH-
%  CENTRED FOR THE MOON) ECLIPTIC SYSTEM.
%  
%  [E, ksun] = uplanet (mjd2000, IBODY);
%
%  These ephemerides were succesfully compared with JPL/NAIF/SPICE
%  ephemerides using de405.bps.
%
%  OUTPUT:
%          E(6) = MEAN KEPLERIAN ELEMENTS OF DATE
%                 E = [a e i Om om theta] [km, rad]
%          ksun = GRAVITY CONSTANT OF THE SUN or OF THE EARTH
%  INPUT :
%          mjd2000 = MODIFIED JULIAN DAY SINCE 01/01/2000, 12:00 (MJD2000)
%                    (MJD2000 = MJD-51544.5)
%          IBODY = NUMBER OF CELESTIAL BODY (< 11)
%                   1:   Mercury
%                   2:   Venus
%                   3:   Earth
%                   4:   Mars
%                   5:   Jupiter
%                   6:   Saturn
%                   7:   Uranus
%                   8:   Neptune
%                   9:   Pluto
%                   10:  Sun
%
% Note: The ephemerides of the Moon are given by EphSS_kep, according to
%       to the algorithm in moon_eph.m
%
%  FUNCTIONS CALLED: none
%
%  - P. Dysli - 1977
% 
%  CHANGELOG
%  - Camilla Colombo - 28/12/06 - tidied up
%  Revised by Matteo Ceriotti - 10/01/2007
%  - Camilla Colombo - 03/05/2008 - case 11 deleted.
%  - Matteo Ceriotti, Camilla Colombo - 11/09/2008 - All ephemerides
%       shifted 0.5 days back in time. Now mjd2000 used in this function is
%       referred to 01/01/2000 12:00. In the old version it was referred to
%       02/01/2000 00:00.
%  - Matteo Ceriotti, Camilla Colombo - 11/09/2008 - Corrected ephemeris of
%       Pluto.
%
%---------------------------------------------------------------------

PI2=2*pi;
RAD=pi/180;

G=6.67259e-20;
msun=1.988919445342813e+030;
ksun=msun*G;

KM  = 149597870.66;

%  T = JULIAN CENTURIES SINCE 31/12/1899 at 12:00
T   = (mjd2000 + 36525)/36525.00;
TT  = T*T;
TTT = T*TT;
E(1)=T*0;

%
%  CLASSICAL PLANETARY ELEMENTS ESTIMATION IN MEAN ECLIPTIC OF DATE
%
switch round(IBODY)
    
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
        E(1) = 39.481686778174627;
        E(2) = 2.4467e-001;
        E(3) = 17.150918639446061;
        E(4) = 110.27718682882954;
        E(5) = 113.77222937912757;
        XM   = 4.5982945101558835e-008;
        E(6) = 1.5021e+001 + XM*mjd2000*86400;
    %
    %  SUN
    %
    case 10
        E = [0 0 0 0 0 0];
    %
    %  MOON (AROUND EARTH)
    %
%     case 11
%         msun=59.736e23;
%         ksun=msun*G;
%         E(1) = 0.0025695549067660;
%         E(2) = 0.0549004890;
%         E(3) = 5.145396388888888890;
%         E(4) = 2.59183275e+2  - 1.93414200833333333e+3*T + 2.07777777777777778e-3*TT + 2.22222222222222222e-6*TTT;
%         E(4) = mod(E(4) + 108e3, 360);
%         E(5) = 7.51462805555555556e+1    + 6.00317604166666667e+3*T - 1.24027777777777778e-2*TT - 1.47222222222222222e-5*TTT;
%         E(5) = mod(E(5), 360);
%         XM   = 4.77198849108333333e+5    + 9.19166666666666667e-3*T  + 1.43888888888888889e-5*TT;
%         E(6) = 2.96104608333333333e2     + XM                    *T;
    
    otherwise
	disp(IBODY)
    if round(IBODY)==11
        error('no planet in the list. For the Moon use EphSS_kep instead')
    else
        error('no planet in the list')
    end
end
  
%
%  CONVERSION OF AU INTO KM, DEG INTO RAD AND DEFINITION OF  XMU
%
 
E(1)   = E(1)*KM;       % a [km]
E(3:6) = E(3:6)*RAD;    % Transform from deg to rad
E(6)   = mod(E(6),PI2);
% XMU  = (XM*RAD/(864000*365250))^2*E(1)^3;
phi    = E(6);          % phi is the eccentric anomaly, uses E(6)=M as a first guess
     
for i=1:5
    g=E(6)-(phi-E(2)*sin(phi)); 
	g_primo=(-1+E(2)*cos(phi));
 	phi=phi-g/g_primo;   % Computes the eccentric anomaly E
end

theta=2*atan(sqrt((1+E(2))/(1-E(2)))*tan(phi/2));

E(6)=theta;

return