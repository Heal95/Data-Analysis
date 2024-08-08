#!/usr/bin/env python

#module Iterate
#coding=utf-8

# Import modules
import Loc_eq_f as lef
import os
import numpy as np

# Model info (Balkan earth model) - parameters are appropriate for Dinarides region
h1 = 30; vp1 = 5.40; vs1 = 3.05
h2 = 16; vp2 = 6.66; vs2 = 3.83
vp3 = 8.00; vs3 = 4.60

def read_bulletin(in_folder, file):

    # Read proposed origin time of the event
    eq_hyp_time = float(file[-4:-2])*3600 + float(file[-2:])*60

    # Read proposed hypocenter location of the event (must be provided in Bulletin file)
    F = open(os.path.join(in_folder,file), 'r')
    lines = F.readlines()
    phi, lamb = lines[1].split() 
    phi_eq_pp = float(phi); lamb_eq_pp = float(lamb)
   
    # Read bulletin files - note: phase names need to be connected with station names
    T = []; sta_phase = [];
    for k in range(2,len(lines)):
        sta_ph, ph_time = lines[k].split()  
        sta_phase.append(sta_ph); T.append(eq_hyp_time+float(ph_time))
    F.close()
    
    return eq_hyp_time, phi_eq_pp, lamb_eq_pp, T, sta_phase

def sta_azimuth_delta(sta_phase, phi_eq_pp, lamb_eq_pp):
    # Define list with stations for which there are phases
    stations = [];
    for k in range(0,len(sta_phase)):
        stations.append(sta_phase[k][:-2])
    
    # Create distance and azimuth lists for used stations 
    delta_sta = []; azimuth_sta = [];
    for sta in stations:
        F = open('./stacoord.lst', 'r')
        lines = F.readlines()
        for line in lines:
            st, phi_sta, lamb_sta, _ = line.split()
            if sta == st:
                D = lef.delta(phi_eq_pp,float(phi_sta),lamb_eq_pp,float(lamb_sta))
                AZ = lef.azimuth(float(phi_sta),phi_eq_pp,lamb_eq_pp,float(lamb_sta),D)
                delta_sta.append(D); azimuth_sta.append(AZ)
        F.close()

    return(delta_sta, azimuth_sta, stations)

def fix_origin_time(sta_phase, T, eq_hyp_time, delta_sta, Depth):
    
    # Fix proposed origin time (calculated from Pg phases of the two closest stations)
    t_pom = []; k_pom = [];
    for k in range(0,len(sta_phase)):
        if sta_phase[k][-2:] == 'PG':
            t_pom.append(T[k]-eq_hyp_time)
            k_pom.append(k)
    tt_pom = t_pom.copy()
    t_pom.sort()
    index1 = tt_pom.index(t_pom[0]); index2 = tt_pom.index(t_pom[1]);
    k1 = k_pom[index1]; k2 = k_pom[index2];
    
    t1 = t_pom[0] - lef.tg(vp1,delta_sta[k1],Depth)
    t2 = t_pom[1] - lef.tg(vp1,delta_sta[k2],Depth)
    eq_hyp_time_pp = eq_hyp_time + (t1+t2)/2

    return(eq_hyp_time_pp)

def Geiger_first_it(sta_phase, T, eq_hyp_time_pp, phi_eq_pp, lamb_eq_pp, delta_sta, azimuth_sta, Depth):
    
    # First iteration of modified Geiger method
    L = []; B = [];
    for k in range(0,len(T)):
        di = 1
        if (sta_phase[k][-2:])=='PG':
            ti = eq_hyp_time_pp + lef.tg(vp1, delta_sta[k], Depth)
            li = ti - T[k]
            ai = lef.aPg(vp1, delta_sta[k], Depth, azimuth_sta[k], phi_eq_pp)
            bi = lef.bPg(vp1, delta_sta[k], Depth, azimuth_sta[k])
            ci = lef.cPg(vp1, delta_sta[k], Depth)
        elif (sta_phase[k][-2:])=='SG':
            ti = eq_hyp_time_pp + lef.tg(vs1, delta_sta[k], Depth)
            li = ti - T[k]
            ai = lef.aPg(vs1, delta_sta[k], Depth, azimuth_sta[k], phi_eq_pp)
            bi = lef.bPg(vs1, delta_sta[k], Depth, azimuth_sta[k])
            ci = lef.cPg(vs1, delta_sta[k], Depth)
        elif (sta_phase[k][-2:])=='PN':
            ti = eq_hyp_time_pp + lef.tn(vp1, vp2, vp3, h1, h2, delta_sta[k], Depth)
            li = ti - T[k]
            ai = lef.aPn(vp3, azimuth_sta[k], phi_eq_pp)
            bi = lef.bPn(vp3, azimuth_sta[k])
            ci = lef.cPn(vp2, vp3)
        else:
            ti = eq_hyp_time_pp + lef.tn(vs1, vs2, vs3, h1, h2, delta_sta[k], Depth)
            li = ti - T[k]
            ai = lef.aPn(vs3, azimuth_sta[k], phi_eq_pp)
            bi = lef.bPn(vs3, azimuth_sta[k])
            ci = lef.cPn(vs2, vs3)
    
        L.append([li])
        B.append([float(ai), float(bi), float(ci), float(di)])
        
    L = np.asarray(L)
    B = np.asarray(B); Bt = B.transpose()
    Q = np.linalg.inv(Bt.dot(B))
    M = Bt.dot(L)
    X = - Q.dot(M)
    
    # Fix values after first iteration
    lamb_eq_pp += X[0]/111.2
    phi_eq_pp += X[1]/111.2
    Depth += X[2]
    eq_hyp_time_pp += X[3]
    i = 0

    return(X, lamb_eq_pp, phi_eq_pp, Depth, eq_hyp_time_pp, i)

def Geiger_it(X, sigma, i, stations, phi_eq_pp, lamb_eq_pp, sta_phase, T, eq_hyp_time_pp, Depth):
    # Iterate loop until desired accuracy is achieved
    while (np.absolute(X[0])>sigma) or (np.absolute(X[1])>sigma) or (np.absolute(X[2])>sigma) or (np.absolute(X[3])>sigma):
        delta_sta = []; azimuth_sta = [];
        for sta in stations:
            F = open('./stacoord.lst', 'r')
            lines = F.readlines()
            for line in lines:
                st, phi_sta, lamb_sta, _ = line.split()
                if sta == st:
                    D = lef.delta(phi_eq_pp,float(phi_sta),lamb_eq_pp,float(lamb_sta))
                    AZ = lef.azimuth(float(phi_sta),phi_eq_pp,lamb_eq_pp,float(lamb_sta),D)
                    delta_sta.append(D); azimuth_sta.append(AZ)
            F.close()
        
        L = []; B = [];
        for k in range(0,len(T)):
            di = 1
            if (sta_phase[k][-2:])=='PG':
                ti = eq_hyp_time_pp + lef.tg(vp1, delta_sta[k], Depth)
                li = ti - T[k]
                ai = lef.aPg(vp1, delta_sta[k], Depth, azimuth_sta[k], phi_eq_pp)
                bi = lef.bPg(vp1, delta_sta[k], Depth, azimuth_sta[k])
                ci = lef.cPg(vp1, delta_sta[k], Depth)
            elif (sta_phase[k][-2:])=='SG':
                ti = eq_hyp_time_pp + lef.tg(vs1, delta_sta[k], Depth)
                li = ti - T[k]
                ai = lef.aPg(vs1, delta_sta[k], Depth, azimuth_sta[k], phi_eq_pp)
                bi = lef.bPg(vs1, delta_sta[k], Depth, azimuth_sta[k])
                ci = lef.cPg(vs1, delta_sta[k], Depth)
            elif (sta_phase[k][-2:])=='PN':
                ti = eq_hyp_time_pp + lef.tn(vp1, vp2, vp3, h1, h2, delta_sta[k], Depth)
                li = ti - T[k]
                ai = lef.aPn(vp3, azimuth_sta[k], phi_eq_pp)
                bi = lef.bPn(vp3, azimuth_sta[k])
                ci = lef.cPn(vp2, vp3)
            else:
                ti = eq_hyp_time_pp + lef.tn(vs1, vs2, vs3, h1, h2, delta_sta[k], Depth)
                li = ti - T[k]
                ai = lef.aPn(vs3, azimuth_sta[k], phi_eq_pp)
                bi = lef.bPn(vs3, azimuth_sta[k])
                ci = lef.cPn(vs2, vs3)
    
            L.append(li)
            B.append([float(ai), float(bi), float(ci), float(di)])
            
        L = np.asarray(L)
        B = np.asarray(B); Bt = B.transpose()
        Q = np.linalg.inv(Bt.dot(B))
        M = Bt.dot(L)
        X = - Q.dot(M)
        
        lamb_eq_pp += X[0]/111.2
        phi_eq_pp += X[1]/111.2
        Depth += X[2]
        eq_hyp_time_pp += X[3]
        i += 1
        if i > 100:
            break
    
    return(eq_hyp_time_pp, phi_eq_pp, lamb_eq_pp, Depth)