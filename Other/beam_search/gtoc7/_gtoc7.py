import numpy as np
from math import *


def asteroid_kep2car(asteroid,mu):
    kep=asteroid[2:8]
    new_asteroid = asteroid

    a   = kep[0];
    e   = kep[1];
    i   = kep[2];
    Om  = kep[4];
    om  = kep[3];
    M   = kep[5];
    tho = mtotheta(M,e);

    rotmat = np.zeros((3,3));
    
    rotmat[0][0] = cos(om)*cos(Om)-sin(om)*cos(i)*sin(Om);
    rotmat[1][0] = cos(om)*sin(Om)+sin(om)*cos(i)*cos(Om);
    rotmat[2][0] = sin(om)*sin(i);

    rotmat[0][1] = -sin(om)*cos(Om)-cos(om)*cos(i)*sin(Om);
    rotmat[1][1] = -sin(om)*sin(Om)+cos(om)*cos(i)*cos(Om);
    rotmat[2][1] = cos(om)*sin(i);

    rotmat[0][2] = sin(i)*sin(Om);
    rotmat[1][2] = -sin(i)*cos(Om);
    rotmat[2][2] = cos(i);

    if ((e<(1+1e-10)) & (e>(1-1e-10))):
        p = 2*a;
    else:
        p = a*(1.0-pow(e,2.0));

    r       = p/(1.0+e*cos(tho));
    xp      = r*cos(tho);
    yp      = r*sin(tho);
    wom_dot = sqrt(mu*p)/pow(r,2.0);
    r_dot   = sqrt(mu/p)*e*sin(tho);
    vxp     = r_dot*cos(tho)-r*sin(tho)*wom_dot;
    vyp     = r_dot*sin(tho)+r*cos(tho)*wom_dot;

    new_asteroid[2] = rotmat[0][0]*xp + rotmat[0][1]*yp;
    new_asteroid[3] = rotmat[1][0]*xp + rotmat[1][1]*yp;
    new_asteroid[4] = rotmat[2][0]*xp + rotmat[2][1]*yp;

    new_asteroid[5] = rotmat[0][0]*vxp + rotmat[0][1]*vyp;
    new_asteroid[6] = rotmat[1][0]*vxp + rotmat[1][1]*vyp;
    new_asteroid[7] = rotmat[2][0]*vxp + rotmat[2][1]*vyp;

    return new_asteroid;

def mtotheta(M,e):
    err = 1.0
    toll = 1.e-7
    E = 0.0
    ratio = 0.0
    if (M < pi):
        E = M + e/2
    else:
        E = M - e/2
    
    while (err > toll):
        ratio = (M - E + e*sin(E) ) / (1 - e*cos(E))
        E = E + ratio
        err = abs(ratio)
    
    tan_theta2 = sqrt((1+e) / (1-e))*tan(E/2)
    theta = 2*atan(tan_theta2)
     
    while (theta < 0  ):
        theta = theta + 2*pi
    
    return theta