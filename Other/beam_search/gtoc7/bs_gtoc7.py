from __future__ import print_function, division
# from _gtoc7 import *
from math import *
from PyKEP import *
import numpy as np
from scipy.spatial import cKDTree as kdtree
from scipy.optimize import fmin_slsqp
import cPickle as pickle # in python 3 no cPickle use pickle
import os
#!/usr/env python

###############################################################################################
###############################################################################################
### beam_search application to finding daughter sequences in gtoc 7 problem. use carefully ;)
### author : carlos ortega absil (carlos.ortega@strath.ac.uk)
###############################################################################################
###############################################################################################




class environment:                          # represents the problem #,
    def __init__(self, ocp = True,
        initial_time_jd=62000.0, mission_time_jd = 6.0*365.25, 
        mass_max_kg = 2000.0, mass_min_kg = 800.0,
        Isp = 3000.0, max_thrust_N=0.3, waiting_time_jd=30.0,
        ephemerides_resolution_jd = 8.0, dt_orbital_distance_scaling_jd = 90.0,
        branching_factor=20, beaming_factor=20, verbosity=1):
        
        compute_or_load_eph = True
        save_load_eph = False # try to save the ephemerides and load the ephemerides if previously saved for exact same initial time and resolution

        self.best = []                          # it stores the best
        self.branching_factor = branching_factor
        self.beaming_factor = beaming_factor

        self.sequence_repr = [962, 1148, 9212, 16075, 3016, 2993, 2341, 6511, 13277, 8442, 5039, 14176]

        # # to get arguments to reproduce sequence stored in file, uncomment lines in select_candidates 
        # # and in pickle.dump(x_sol...) and in main set a0 to be the first of the sequence
        # # also disable .load ephemerides for speed

        # FIXED ATTRIBUTES OF THE ENVIRONMENT  
        # doublecheck all the nondimensionalisations      
        # abuse of the subscript _si sometimes is km 
        self.ocp = ocp

        self.mu = 1.0
        self.mu_si = 1.32712440018e11
        self.au_si = 1.49597870691e8
        self.g0_si = 9.80665
        self.tu_si = sqrt(self.au_si/self.mu_si)*self.au_si
        self.secondsinaday = 86400.0
        self.days2du = sqrt(self.mu_si/self.au_si)*self.secondsinaday/self.au_si
        self.du2days = 1.0/self.days2du
        print('AU:: ',self.au_si)
        print('TU:: ',self.tu_si)
        print('VU:: ',self.au_si/self.tu_si)

        self.tw = waiting_time_jd*self.days2du
        self.initial_time = initial_time_jd*self.days2du
        self.final_time = (initial_time_jd+mission_time_jd)*self.days2du
        self.dt_orbital_distance_scaling = dt_orbital_distance_scaling_jd*self.days2du # This parameter feels important, there should b a log of how it relates to am etc.
        print('OD scaling:: ',self.dt_orbital_distance_scaling)
        self.mass_max_kg = mass_max_kg
        self.mass_min_kg = mass_min_kg
        self.Ispg0 = Isp*self.g0_si*self.secondsinaday/self.days2du/self.au_si/1e3
        self.Tm = max_thrust_N*self.secondsinaday*self.secondsinaday/self.days2du/self.days2du/self.au_si/1e3
        self.Isp_si = Isp
        self.Tm_si = max_thrust_N
        # print('am:::::::::',self.am)
        # read asteroids
        self.listofasteroids = self._read_asteroids()
        self.n_asteroids = len(self.listofasteroids) # includes the earth
        self.ephemerides_resolution = ephemerides_resolution_jd * self.days2du

        if compute_or_load_eph:
            # compute the ephemerides of all asteroids as well as KDTrees for discrete times and store
            self.ephemerides_trees = []
            self.ephemerides_lists = []        
            eph_filename = "ephemerides_initial_"+str(initial_time_jd)+"_final_"+str(initial_time_jd+mission_time_jd)+"_resolution_"+str(ephemerides_resolution_jd)+"_scaling_"+str(dt_orbital_distance_scaling_jd)+".p";
            if save_load_eph and os.path.isfile(eph_filename):
                print('start loading ephemerides')
                eph_file = open( eph_filename, "rb" )
                # self.ephemerides_lists= pickle.load(eph_file)
                print('done loading ephemerides')
            else:
                t = self.initial_time
                print('start precomputing ephemerides')
                while t<=(self.final_time+max(2*self.tw,2*self.ephemerides_resolution)): # factors for safety not to acces uninitialized trees
                    ephemerides_list = self.get_asteroid_ephemerides('all',t,dt_orbital_distance_scaling = self.dt_orbital_distance_scaling)
                    # print('-'*20+'\n',np.asarray(ephemerides_list[:10]))
                    self.ephemerides_lists.append(ephemerides_list)
                    # self.ephemerides_times.append(t)
                    t += self.ephemerides_resolution
                print('done precomputing ephemerides')
                if save_load_eph:
                    eph_file = open(eph_filename, "wb" )
                    pickle.dump(self.ephemerides_lists, eph_file )
                    print('done saving ephemerides')
            print('start putting ephemerides in KDTree form')
            for all_asteroids_eph in self.ephemerides_lists:
                self.ephemerides_trees.append(kdtree(all_asteroids_eph,leafsize = 16))
            print('done putting ephemerides in KDTree form')

        self.verbosity =   verbosity

        self.depth = 0                          # and attributes that change according to depth

    def _read_asteroids(self):
        fileread = open("gtoc7_asteroids_data.txt",'r').readlines()
        listofasteroids = []
        for line in fileread[2:]:
            asteroid_data = line.split()
            asteroid = [int(asteroid_data[0])]+[float(asteroid_data[1])*self.days2du]+[float(x) for x in asteroid_data[2:4]]+[float(x)*pi/180 for x in asteroid_data[4:8]]+[" ".join(asteroid_data[8:])]
            listofasteroids.append(asteroid)
        listofasteroids_car = [asteroid_kep2car(asteroid,1.0) for asteroid in listofasteroids]
        # print('some asteroids :\n',np.array(listofasteroids_car[0:-1:1000]))
        return listofasteroids_car
    
    def get_asteroid_ephemerides(self,asteroid_id,t,dt_orbital_distance_scaling=0.0): # takes t in TU
        # the dt for orbital scaling is a dt that estimates the time of transfer in the linearised lambert problem
        # if dt == 0 then the the adimensionalized state vector [r,v] with AU and AU/TU is returned
        squeeze_result = False
        if asteroid_id == 'all':
            listofasteroids = self.listofasteroids
            numberofasteroids = len(self.listofasteroids)
        elif isinstance(asteroid_id,list):
            listofasteroids = [self.listofasteroids[i] for i in asteroid_id]
            numberofasteroids = len(asteroid_id)
        elif isinstance(asteroid_id,int):
            listofasteroids = [self.listofasteroids[asteroid_id]]
            numberofasteroids = 1
            squeeze_result = True
        else:
            print('ERROR: trying to get ephemerides with invalid asteroid_id', asteroid_id)
        # ephemerides = [[0,0,0,0,0,0]]*numberofasteroids
        ephemerides=[]
        count = 0
        # print(asteroid_ids)
        for asteroid in listofasteroids:
            rf, vf = propagate_lagrangian(asteroid[2:5],asteroid[5:8],t-asteroid[1],1.0)
            if dt_orbital_distance_scaling>1e-6:
                rdt = np.array(rf)/(dt_orbital_distance_scaling)
                rf = rdt + np.array(vf)
                vf = rdt
            ephemerides.append(list(rf)[:]+list(vf)[:])
            # ephemerides[count][:3] = list(rf)[:3]
            # ephemerides[count][3:] = list(vf)[3:]
            # count+=1
            # if count < 10 : print('*'*20+'\n',np.asarray(ephemerides[count]))
        # print('-'*20+'\n',np.asarray(ephemerides[:10]))
        return ephemerides[0] if squeeze_result else ephemerides

    def update_best(self,level):            # this method updates the best
        # here the idea would be to take a list with the non.dominated in rm and rt
        # for the moment just the minimum rt
        if level:
            self.best = max(level,key=lambda x:x.rt)

    def update_environment(self,depth):     # this method updates the environment according to depth
        self.depth = depth


    def beam(self,level):
        if self.beaming_factor>0:
            lt = sorted(level,key= lambda x:-x.rt)[:self.beaming_factor]
            # # mix mass and time
            # 
            # lm = sorted(level,key= lambda x:-x.rm)[:self.beaming_factor]
            # nextlevel = []
            # i=0
            # while len(nextlevel)<=self.beaming_factor and i<min(len(lt),len(lm)):
            #     r = np.random.rand()
            #     if r<0.85:
            #         if lt[i] not in nextlevel:
            #             nextlevel.append(lt[i])
            #     else:
            #         if lm[i] not in nextlevel:
            #             nextlevel.append(lm[i])
            #     i+=1
            # return nextlevel[:self.beaming_factor] 
            # set the beaming factor for time as recommended
            return lt 
        else:
            return level

    def beam_search(self,level,max_depth=10000):
        if isinstance(level,list):
            b = level # b is this generation, level is the initial set or initial node
        else:
            print('ERROR: need to pass a list of nodes to initialise beam_search')
        depth = 0

        self.update_best(level)                  # update best sequence
        print('bestguy',self.best.rt)
        while len(b)>0 and depth<max_depth:                          # while there is a next generation
            depth+=1
            l=[]
            for n in b:                          # take te nodes
                l+=branch(n, self)               # make them bear children
            
            self.update_best(l)                  # update the environment
            self.update_environment(depth)       

            if self.verbosity>0: print('size of level before beaming : ', len(l))

            b = self.beam(l)                     # select some beams, here is where the termination condition dwindles

            if self.verbosity>0: print('size of level after beaming : ', len(b))
            if self.verbosity>0: print('BEST: length : ',len(self.best.sequence),'  / t0 : ', self.best.t0*self.du2days,'  / m0 : ', self.best.m0 , '\n sequence :', self.best.sequence)
            pickle.dump(sorted(l,key= lambda x:-x.rt), open('last_level962.p','wb'))


class node:
    def __init__(self,sequence,t0,m0,environment):      # t0 in DU, m0 in kg
        self.sequence = sequence
        # store the basic attributes passed as arguments
        self.asteroid_id = sequence[-1]
        self.t0 = t0
        self.m0 = m0
        # compute the remaining budgets based on environment
        self.rt = (environment.final_time-t0)/(environment.final_time-environment.initial_time)
        self.rm = (m0-environment.mass_min_kg)/(environment.mass_max_kg-environment.mass_min_kg)

    def give_birth(self,environment):
        if self.sterile(environment): return []
        children_id = self.select_candidates(environment)
        children = []
        for kid_id in children_id:
            if environment.ocp:
                child_args = self.compute_child_args_sf(kid_id,environment) # in this case the costly computations are here
            else:
                child_args = self.compute_child_args_lambert(kid_id,environment) # in this case the costly computations are here
            # note that several children should be born not only for minimum time!
            if self.is_child_legal(child_args,environment): children.append(child_args)
        return children

    def sterile(self,environment):
        return (self.t0+environment.tw>=environment.final_time) or (self.m0<environment.mass_min_kg + 10.0)# the node is sterile if it cant wait at next asteroid before end of mission

    def select_candidates(self,environment):         # returns a list of next level asteroids possible to visit as in the gtoc paper
        # return [environment.sequence_repr[environment.depth+1]]
        i_t_tree = int(round((self.t0-environment.initial_time)/environment.ephemerides_resolution))  # put t0 to tree resolution
        # NOTE you are querying the time at t0 so dismissing waiting time
        current_eph = environment.ephemerides_lists[i_t_tree][self.asteroid_id]             # get the ephemerides of you in i_t
        eph_tree = environment.ephemerides_trees[i_t_tree]                                  # get the tree to query
        # query the tree for the nearest neighbors with a branching factor:
        distances, nearby_id = eph_tree.query(current_eph,k=environment.branching_factor+len(self.sequence)+2)
        children = [x for x in nearby_id if (x not in self.sequence) and x>0 and x<environment.n_asteroids]
        return children[:environment.branching_factor] 
        # Note that select candidates only clusters for t0

    def is_child_legal(self,child_args,environment):# check if the child is legal before initializing it
        # is legal if the budgets are positive simply, necessary to avoid dead kids
        return child_args and child_args[1]<=environment.final_time and child_args[2]>=environment.mass_min_kg

    def compute_child_args_lambert(self,kid_id,environment):    # do the costly computations
        # solve the optimal control problem here or pause and ask Gianluca to solve it by hand
        # start by applying the OCP logic to the lambert with 2 delta V and test the global before going on...

        # with lambert ONLY CONSIDERING MINIMUM TIME TRANSFERS otherwise put more 
        tof_max = 365.25*environment.days2du # upper bound on tf-t0
        func   = lambda x:x[1] # t_f to minimise
        cfun  = lambda x : self.cfun_lambert(kid_id,x[0],x[1],environment) # dv and tof constraints
        bounds = [(self.t0, self.t0+tof_max),(self.t0, self.t0+tof_max)] # bounded btw t0 and t0 + some days
        x0 = np.array([self.t0 , self.t0+1.0*environment.days2du])
        ieqcons = [lambda x: x[1] - x[0] -1e-6, cfun]
        x= fmin_slsqp(func, x0, bounds = bounds, ieqcons=ieqcons,full_output =False, iprint =environment.verbosity-1)
        # print('opt_sol = ', x)
        # print('ieqcons = ', [f(x) for f in ieqcons])

        child_args = []
        # recheck constraints for safety
        if x[1]>x[0]:
            dv = self.dv_lambert(kid_id,x[0],x[1],environment)
            if self.cfun_lambert(kid_id,x[0],x[1],environment)>0:
                # print('dv = ',dv)
                child_sequence = self.sequence + [kid_id]
                child_t0 = x[1] + environment.tw
                # print(self.m0,dv,environment.Ispg0)
                child_m0 = self.m0/exp(dv/environment.Ispg0) # rocket equation
                child_args=[child_sequence, child_t0, child_m0]
                if environment.verbosity>1:
                    print('successful transfer! : ', self.t0*environment.du2days, x[0]*environment.du2days, x[1]*environment.du2days, dv*environment.au_si/environment.tu_si)
                    print('a0 :',environment.get_asteroid_ephemerides(self.asteroid_id,x[0],0.0)[:])
                    print('a1 :',environment.get_asteroid_ephemerides(kid_id,x[1],0.0)[:])
        return child_args

    def compute_child_args_sf(self,kid_id,environment):    # do the costly computations
        # solve the optimal control problem here or pause and ask Gianluca to solve it by hand
        # start by applying the OCP logic to the lambert with 2 delta V and test the global before going on...

        # ONLY CONSIDERING MINIMUM TIME TRANSFERS otherwise put more 
        n_seg = 10
        tofmax = 365.0

        t0 = self.t0*environment.du2days
        # x will be [ts, tf, mf, Tx1, Ty1, Tz1, Tx2.... Tzn]. length of x defines number of segments, ts and tf in day units
        
        # func   = lambda x:x[1] # t_f to minimise
        # func   = lambda x:x[2] # m_f to minimise
        func   = lambda x:-(0.5*(environment.final_time-x[1]*environment.days2du)/(environment.final_time-environment.initial_time)+0.5*(x[2]-environment.mass_min_kg)/(environment.mass_max_kg-environment.mass_min_kg))
        
        x0 = [t0 , t0 + 60.0 , max(self.m0 - 50.0, environment.mass_min_kg+1.0)]

        # # dumb initial guess for T... full throttle in the direction of the asteroid you target
        # ephs = environment.get_asteroid_ephemerides(self.asteroid_id,self.t0,0.0)
        # ephf = environment.get_asteroid_ephemerides(kid_id,self.t0+60.0*environment.days2du,0.0)
        # rs = np.array(ephs[:3])
        # rf = np.array(ephf[:3])

        # dr = rf-rs
        # direction = 0.99*dr/np.linalg.norm(dr)

        ## initial guess with lambert and a time of flight of 100 days 
        (dv1,dv2) = self.dv2_lambert(kid_id,self.t0,self.t0+60.0*environment.days2du,environment)
        # divide segments in dv1 and dv2
        dv1norm = np.linalg.norm(dv1)
        dv2norm = np.linalg.norm(dv2)
        thrust1 = np.array(dv1)/dv1norm
        thrust2 = np.array(dv2)/dv2norm
        n_seg1 = round( n_seg * dv1norm / (dv1norm + dv2norm))
        x0.extend(list(thrust1)* (int (n_seg1)) + list(thrust2)*(int (n_seg-n_seg1)))
        # # all segments equal
        # thrust = np.array(dv1)+np.array(dv2)
        # thrust /= np.linalg.norm(thrust)
        # x0.extend(list(thrust)*n_seg)
        # # print(x0)



        ceqfun = lambda x: self.ceqfun_sims_flanagan(kid_id,x,environment)
        cineqfun  = lambda x : self.cineqfun_sims_flanagan(kid_id,x,environment) # dv and tof constraints
        
        bounds = [(t0, t0+tofmax-0.5),(t0+0.5, t0+tofmax),(max(self.m0-300.0, environment.mass_min_kg+0.5),self.m0)] # bounded btw t0 and t0 + some days
        bounds.extend([(-2.0,2.0)]*3*n_seg)


        # last thing I try forget optimality just minimize ineq. constraint square sum
        func_eq = lambda x: sum([i*i for i in ceqfun(x)])
        x02 = fmin_slsqp(func_eq,x0,bounds=bounds,f_ieqcons=cineqfun,full_output=False,iprint=environment.verbosity-1,acc=1e-5,iter=500)
        x_sol = x02[:]
        if func_eq(x02) < 5.0e-2:
            x_sol= fmin_slsqp(func, x02, bounds = bounds,  f_ieqcons=cineqfun, f_eqcons=ceqfun, full_output =False, iprint =environment.verbosity-1,acc=1e-5,iter=500)
        # print('opt_sol = ', x_sol)

        child_args = []

        child_sequence = self.sequence[:] + [kid_id]
        child_t0 = x_sol[1]*environment.days2du + environment.tw
        child_m0 = x_sol[2]
        # print(x_sol)
        # pickle.dump(x_sol,open('x'+str(environment.depth+1)+'.p','wb'))
        # recheck constraints for safety
        if all([abs(i)<1e-5 for i in ceqfun(x_sol)]) and all([i>=-1e-5 for i in cineqfun(x_sol)]) and all([(bound[0]-1e-6 <= xi <= bound[1]+1e-6) for (xi,bound) in zip(x_sol,bounds)]):
            child_args=[child_sequence, child_t0, child_m0]
            if environment.verbosity>1:
                print('successful transfer! ::: t0 = ', t0, ' /  ts = ', x_sol[0] , ' /  tf = ', x_sol[1], ' /  m0 = ', self.m0, ' /  mf = ', x_sol[2])
                print('a0 :',self.asteroid_id,' / ',  environment.get_asteroid_ephemerides(self.asteroid_id,x_sol[0]*environment.days2du,0.0)[:])
                print('a1 :',kid_id, ' / ', environment.get_asteroid_ephemerides(kid_id,x_sol[1]*environment.days2du,0.0)[:])
                print('control = ',x_sol[3:])
                print('ceq(xF)   =',ceqfun(x_sol))
                print('cineq(xF) =',cineqfun(x_sol))
        elif environment.verbosity>1:
            print('ABORTING KID, UNSUCCESSFUL OPT\n','======================\n','t0 = ', t0, ' /  ts = ', x_sol[0] , ' /  tf = ', x_sol[1], ' /  m0 = ', self.m0, ' /  mf = ', x_sol[2])
            print('ceq(x0)   =',ceqfun(x0))
            print('cineq(x0) =',cineqfun(x0))
            print('ceq(xF)   =',ceqfun(x_sol))
            print('cineq(xF) =',cineqfun(x_sol))
            print('a0 :',self.asteroid_id,' / ',  environment.get_asteroid_ephemerides(self.asteroid_id,x_sol[0]*environment.days2du,0.0)[:])
            print('a1 :',kid_id, ' / ', environment.get_asteroid_ephemerides(kid_id,x_sol[1]*environment.days2du,0.0)[:])
        return child_args

    def ceqfun_sims_flanagan(self,a1_id,x,environment):
        leg = self.compute_leg_sims_flanagan(a1_id,x,environment)
        mm = leg.mismatch_constraints()
        res = list(np.array(mm[:3])/environment.au_si/1e3)+ list(np.array(mm[3:6])*environment.tu_si/environment.au_si/1e3)+ [mm[6]/self.m0]
        return res

    def cineqfun_sims_flanagan(self,a1_id,x,environment):
        leg = self.compute_leg_sims_flanagan(a1_id,x,environment)
        cineq =  [-i for i in leg.throttles_constraints()]
        cineq.append((x[1]-x[0]-1e-6)/365.25)
        return cineq

    def compute_leg_sims_flanagan(self,a1_id,x,environment):
        ts = x[0]
        tf = x[1]
        mf = x[2]
        throttles = x[3:]
        TS = epoch(ts,'mjd')
        TF = epoch(tf,'mjd')
        ephs = environment.get_asteroid_ephemerides(self.asteroid_id,ts*environment.days2du,0.0)
        ephf = environment.get_asteroid_ephemerides(a1_id,tf*environment.days2du,0.0)
        rs = np.array(ephs[:3])*environment.au_si*1e3
        vs = np.array(ephs[3:])*environment.au_si/environment.tu_si*1e3
        rf = np.array(ephf[:3])*environment.au_si*1e3
        vf = np.array(ephf[3:])*environment.au_si/environment.tu_si*1e3
        sc = sims_flanagan.spacecraft(self.m0,environment.Tm_si,environment.Isp_si)
        xs = sims_flanagan.sc_state(rs,vs,sc.mass)
        xf = sims_flanagan.sc_state(rf,vf,mf)
        # print(rs,vs)
        # print(rf,vf)
        # print(xs)
        # print(xf)
        leg = sims_flanagan.leg(TS,xs,throttles,TF,xf,sc,environment.mu_si*1e9)
        leg.high_fidelity = True
        return leg


    def dv_lambert(self,a1_id,ts,tf,environment):
        ephs = environment.get_asteroid_ephemerides(self.asteroid_id,ts,0.0)
        ephf = environment.get_asteroid_ephemerides(a1_id,tf,0.0)
        rs = ephs[:3]
        vs = ephs[3:]
        rf = ephf[:3]
        vf = ephf[3:]        
        # print(self.asteroid_id,a1_id)
        # print('rs:',len(rs))
        # print('rf:',len(rf))
        # print('t :',tf-ts)
        l  = lambert_problem(rs,rf,tf-ts,1.0,False,0) #note the HARDCODED prograde and multirevolution flags
        v1 = l.get_v1()[0]
        v2 = l.get_v2()[0]
        dv = np.linalg.norm(np.array(v1)-np.array(vs))+np.linalg.norm(np.array(v2)-np.array(vf))
        return dv

    def dv2_lambert(self,a1_id,ts,tf,environment):
        ephs = environment.get_asteroid_ephemerides(self.asteroid_id,ts,0.0)
        ephf = environment.get_asteroid_ephemerides(a1_id,tf,0.0)
        rs = ephs[:3]
        vs = ephs[3:]
        rf = ephf[:3]
        vf = ephf[3:]        
        # print(self.asteroid_id,a1_id)
        # print('rs:',len(rs))
        # print('rf:',len(rf))
        # print('t :',tf-ts)
        l  = lambert_problem(rs,rf,tf-ts,1.0,False,0) #note the HARDCODED prograde and multirevolution flags
        v1 = l.get_v1()[0]
        v2 = l.get_v2()[0]
        dv1 = np.array(v1)-np.array(vs)
        dv2 = np.array(v2)-np.array(vf)
        
        return (dv1,dv2)

    def cfun_lambert(self,a1_id,ts,tf,environment):
        return (environment.Tm*(tf-ts)/self.m0-self.dv_lambert(a1_id,ts,tf,environment)) if tf-ts > 0.0 else -1.0

    def filter_children(self,children,environment):  # this is late abortion at parent level
        return children                              # in this case the children live happily
    

def branch(parent,environment):
    children=[]   
    child_args_list = parent.give_birth(environment)            # returns a list of arguments to initialise new nodes
    for child_args in child_args_list:
        new_node = node(*child_args,environment =environment)   # initialise the node and do the costly computations
        children.append(new_node)                               # appends new node object to the next level
    parent.filter_children(children,environment)                # here is where the branching factor comes into play
    return children

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

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import time

    env = environment(ocp=True)
    t0 = env.initial_time+env.tw
    m0=  2000.0
    # t0_rgs = [(51544.0+7500.0)*env.days2du, (51544.0+12000.0)*env.days2du] # need to play with eph loading to be able to do this
    a0_rgs = [1, 16256]
    # a0 = [10000]

    best = []
    for i in range(1):
        # t0 = t0_rgs[0]+np.random.rand()*(t0_rgs[1]-t0_rgs[0]) # need to play with eph loading to be able to do this
        a0 = [np.random.random_integers(a0_rgs[0],a0_rgs[1])]
        a0 = [962]
        print('a0 : ',a0)
        level0 = [node(a0,t0,m0,env)]
        
        env.beam_search(level0)

        print(env.best.sequence)

        best.append([env.best.sequence,env.best.t0*env.du2days,env.best.m0])

    print ('###############')
    print ('### RESULTS ###')
    print ('===============')
    count = 0
    for i in sorted(best,key=lambda x : -len(x[0])):
        count += 1
        print ('###### ',count,' ######')
        print ('length : ' , len(i[0]) , '  /  tf : ' , i[1] , '  /  mf : ' , i[2] )
        print ('  /  sequence : ' , i[0])
    #### TEST FOR EPHEMERIDES
    # for eph in env.ephemerides_lists:
    #     epha = np.asarray(eph)
    #     ax = plt.figure(1).add_subplot(111, projection='3d')
    #     # ax.scatter(epha[0:10, 0], epha[0:10, 1], epha[0:10, 2])
    #     scat = ax.scatter(epha[:200, 0], epha[:200, 1], epha[:200, 2])



    #     # ax = plt.figure(1).add_subplot(111)
    #     # ax.scatter(eph[0:10, 0], eph[0:10, 1],c = np.transpose(colors[0:10]))


    #     ax.set_xlim(-5, 5)
    #     ax.set_ylim(-5, 5)
    #     ax.set_zlim(-5, 5)
    #     # ax.set_xlabel('x' + str(dimensions[0] + 1))
    #     # ax.set_ylabel('x' + str(dimensions[1] + 1))
    #     # ax.set_zlabel('x' + str(dimensions[2] + 1))
    #     # plt.show(f)
    #     plt.draw()
    #     scat.remove()
    #     # time.sleep(0.5)
    #     # plt.pause(0.05)

    ###### TEST FOR FINDING NEIGHBORS
    if False:
        plt.ion()
        time = 60100.0*env.days2du # DU
        mass = 2000.0 # kg
        unscale = True
        # print(env.days2du)
        # print(env.ephemerides_resolution)
        # print(env.initial_time)
        node_asteroid = 10000
        nod = node([1, 100, 3, 290, node_asteroid], time, mass, env)
        next_asteroids = nod.select_candidates(env)
        print(next_asteroids)
        i_time_eph = int(round((time-env.initial_time)/env.ephemerides_resolution))

        ax = plt.figure(1).add_subplot(111, projection='3d')
        for k in range(100):
            for i in range(i_time_eph-10, i_time_eph+11):
                # get eph
                eph0 = [env.ephemerides_lists[i][j] for j in [node_asteroid]+next_asteroids]
                # unscale
                if unscale and (env.dt_orbital_distance_scaling>1e-6):
                    eph = []
                    ast_eph = [0,0,0,0,0,0]
                    for ast_eph0 in eph0:
                        for i in range(3):
                            ast_eph[3+i] = ast_eph0[i]-ast_eph0[3+i]
                        for i in range(3):
                            ast_eph[i] = ast_eph0[3+i]*env.dt_orbital_distance_scaling
                        eph.append(ast_eph[:])
                else:
                    eph = eph0;
                eph = np.array(eph)
                # ax.scatter(eph[:, 0], eph[:, 1], eph[:, 2])
                me_scat    = ax.scatter(eph[0,0], eph[0,1], eph[0,2],c='red')
                other_scat = ax.scatter(eph[1:, 0], eph[1:, 1], eph[1:, 2],c='black')
                ax.set_xlim(-5, 5)
                ax.set_ylim(-5, 5)
                ax.set_zlim(-5, 5)
                plt.draw()
                plt.pause(0.2)
                me_scat.remove()
                other_scat.remove()
        

    # # set the parameters of the node and initialise node
    # sequence = [(0,0)]
    # fitness = 0
    # initial_node = node(sequence,fitness,env)

    # # do the beam search
    # env.beam_search([initial_node])

    # # print the best
    # print('=======================')
    # print('fitness : ', env.best_fitness)
    # print('n. sol. : ', len(env.best))
    # # for sol in env.best:
    # #     print(sol.visited)
