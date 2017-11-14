from PyKEP import *
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import cPickle as pickle
from bs_gtoc7 import *
from math import *

directory = 'seq_12'
env = environment()
t0 = 62000.0
T0=epoch(t0,'mjd')
m0 = 2000.0

seq=[13414, 1898, 10460, 2654, 8518, 2091, 16002, 16067, 6567, 15984, 406]
seq=[962, 1148, 9212, 16075, 3016, 2993, 2341, 6511, 13277, 8442, 5039, 14176]

fig = plt.figure()
axis = fig.gca(projection='3d')
axis.set_zlim([-.5,.5])

# orbit_plots.plot_planet(planet.gtoc7(seq[0]), T0, units=AU, legend=False, ax = axis)

# T1=epoch(t0+30.0,'mjd')
# T2=epoch(t0+120.0,'mjd')
# orbit_plots.plot_planet(planet.gtoc7(seq[0]), T1, units=AU, legend=False, ax = axis)
# orbit_plots.plot_planet(planet.gtoc7(seq[0]), T2, units=AU, legend=False, ax = axis)

ephs = env.get_asteroid_ephemerides(seq[0],t0*env.days2du,0.0)
ephe = env.get_asteroid_ephemerides(0,t0*env.days2du,0.0)

# ephf = env.get_asteroid_ephemerides(seq[0],(t0+30)*env.days2du,0.0)
axis.scatter([ephs[0]],[ephs[1]],[ephs[2]],c='m',s=40)
# axis.scatter([ephf[0]],[ephf[1]],[ephf[2]],c='r',s=40)

# (xx,yy,zz) = ([],[],[])
# for l in range(100):
#     k = -4.0+ l*8.0/99.0
#     xx.append(k)
#     yy.append(k)
#     zz.append(0.0)

# axis.plot(xx,yy,zz)
# axis.scatter([ephff[0]],[ephff[1]],[ephff[2]],c='r',s=40)

i=0
ephs = env.get_asteroid_ephemerides(seq[0],t0*env.days2du,0.0)

axis.scatter([0],[0],[0],c='y',s=50)
tx = t0*env.days2du
dt = 10.0*env.days2du
xab = [ephs[0]]
yab = [ephs[1]]
zab = [ephs[2]]
xe = [ephe[0]]
ye = [ephe[1]]
ze = [ephe[2]]
for j in range(200):
    tx += dt
    ephx = env.get_asteroid_ephemerides(seq[0],tx,0.0)
    xab.append(ephx[0])
    yab.append(ephx[1])
    zab.append(ephx[2])

    ephex = env.get_asteroid_ephemerides(0,tx,0.0)
    xe.append(ephex[0])
    ye.append(ephex[1])
    ze.append(ephex[2])

axis.plot(xab,yab,zab,c='c')
axis.plot(xe,ye,ze,c='b')
axis.scatter([ephe[0]],[ephe[1]],[ephe[2]],c='b',s=40)

for i in range(1,len(seq)):
    if i==1:
        sc = sims_flanagan.spacecraft(m0,env.Tm_si,env.Isp_si)
    else:
        x_file = open( directory+'/x'+str(i-1)+'.p', "rb" )
        x= pickle.load(x_file)
        sc = sims_flanagan.spacecraft(x[2],env.Tm_si,env.Isp_si)

    x_file = open( directory+'/x'+str(i)+'.p', "rb" )
    x= pickle.load(x_file)
    ts = x[0]
    tf = x[1]
    mf = x[2]
    throttles = x[3:]
    TS = epoch(ts,'mjd')
    TF = epoch(tf,'mjd')


    ephs = env.get_asteroid_ephemerides(seq[i-1],ts*env.days2du,0.0)
    ephf = env.get_asteroid_ephemerides(seq[i],tf*env.days2du,0.0)
    axis.scatter([ephs[0]],[ephs[1]],[ephs[2]],c='g',s=40)
    axis.scatter([ephf[0]],[ephf[1]],[ephf[2]],c='r',s=40)

    tx = tf*env.days2du
    dt = 10.0*env.days2du
    xab = [ephf[0]]
    yab = [ephf[1]]
    zab = [ephf[2]]
    for j in range(200):
        tx += dt
        ephx = env.get_asteroid_ephemerides(seq[i],tx,0.0)
        xab.append(ephx[0])
        yab.append(ephx[1])
        zab.append(ephx[2])

    axis.plot(xab,yab,zab,c='c')

    rs = np.array(ephs[:3])*env.au_si*1e3
    vs = np.array(ephs[3:])*env.au_si/env.tu_si*1e3
    rf = np.array(ephf[:3])*env.au_si*1e3
    vf = np.array(ephf[3:])*env.au_si/env.tu_si*1e3
    
    xs = sims_flanagan.sc_state(rs,vs,sc.mass)
    xf = sims_flanagan.sc_state(rf,vf,mf)

    leg = sims_flanagan.leg(TS,xs,throttles,TF,xf,sc,env.mu_si*1e9)
    leg.high_fidelity = True
    

    orbit_plots.plot_sf_leg(leg, units=AU, N=6, ax=axis)

    # orbit_plots.plot_planet(planet.gtoc7(seq[i-1]), TS, units=AU, legend=False, ax = axis)
    # orbit_plots.plot_planet(planet.gtoc7(seq[i]), TF, units=AU, legend=False, ax = axis)
    axis.scatter([ephs[0]],[ephs[1]],[ephs[2]],c='g',s=40)
    axis.scatter([ephf[0]],[ephf[1]],[ephf[2]],c='r',s=40)
# axis.get_legend().set_visible(False)
plt.show()