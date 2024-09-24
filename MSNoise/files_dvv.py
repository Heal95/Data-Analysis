#!/usr/bin/env python
# coding: utf-8

# Reads dt/t output files from MSNoise and formats it to
# dv/v files for ZZ, RR and TT components

# Import modules
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import os
from datetime import timedelta, date
import pylab
from obspy.geodetics.base import gps2dist_azimuth
import math

# Time axis handeling
def daterange(date1, date2):
    for n in range(int ((date2 - date1).days)+1):
        yield date1 + timedelta(n)

# Median absolute devation
def mad(arr):
    """ Robust version of standard deviation.
        Indices variabililty of the sample.
    """
    arr = np.ma.array(arr).compressed()
    med = np.median(arr)
    return np.median(np.abs(arr - med))

# Read station list
file = open('./Stations.txt',"r")
lines = file.readlines()
file.close()

sta_pairs = []; azs = [];

# Form station pairs and define their azimuths
for k in range(0,len(lines)):
    for j in range(k, len(lines)):
        lat1, lon1, sta1, net1 = lines[k].split()
        lat2, lon2, sta2, net2 = lines[j].split()
        sta_pairs.append(net1+'_'+sta1+'_'+net2+'_'+sta2)
        
        _, az, _ = gps2dist_azimuth(float(lat1),float(lon1),float(lat2),float(lon2))
        azimuth = az * np.pi / 180
        azs.append(azimuth)

# List of data to be considered
# Stacks
days = ['30', '60', '90']
# Start date
start_dt = date(2009, 1, 1)
# End date
end_dt = date(2019, 1, 1)
# Filters
filters = ['0']

# Read original output data and combine it to single file
# ZZ components stays same, RR and TT are calculated from NN, EE, NE and EN
for filt in filters:
    for x in range(0,len(days)):
        # Read dv/v data for ZZ component
        ZZ_sigma, ZZ_med, ZZ_t, ZZ_pts = = [], [], [], []
        j = 0

        for dt in daterange(start_dt, end_dt):
            assist = [];
            exists = os.path.isfile('./DTT/0'+filt+'/0'+days[x]+'_DAYS/ZZ/'+dt.strftime('%Y-%m-%d')+'.txt')
            ZZ_t.append(dt.strftime('%Y-%m-%d')); ZZ_sigma.append(0); ZZ_med.append(0)

            if exists:
                file = open('./DTT/0'+filt+'/0'+days[x]+'_DAYS/ZZ/'+dt.strftime('%Y-%m-%d')+'.txt',"r")
                line = file.readlines()
                file.close()

                for k in range(1,len(line)):
                    Date,Pairs,M,EM,A,EA,M0,EM0 = line[k].split(",")
                    for pair in sta_pairs:
                        if Pairs.strip() == pair:
                            assist.append(float(M))

            ZZ_pts.append(len(assist)); ZZ_sigma[j] = -100*(mad(assist)); ZZ_med[j] = -100*np.median(assist)
            j += 1

        ZZ_sigma, ZZ_med, ZZ_pts = map(np.asarray, [ZZ_sigma, ZZ_med, ZZ_pts])
        
        # Read dv/v to form RR and TT components (azimuthal rotation of NN, EE, EN and NE components)
        RR_sigma, RR_med, RR_t, RR_pts = [], [], [], []
        TT_sigma, TT_med, TT_t, TT_pts = [], [], [], []
        j = 0

        for dt in daterange(start_dt, end_dt):
            assistR = []; assistT = [];
            RR_t.append(dt.strftime('%Y-%m-%d')); RR_sigma.append(0); RR_med.append(0)
            TT_t.append(dt.strftime('%Y-%m-%d')); TT_sigma.append(0); TT_med.append(0)

            ex1 = os.path.isfile('./DTT/0'+filt+'/0'+days[x]+'_DAYS/NN/'+dt.strftime('%Y-%m-%d')+'.txt')
            ex2 = os.path.isfile('./DTT/0'+filt+'/0'+days[x]+'_DAYS/EN/'+dt.strftime('%Y-%m-%d')+'.txt')
            ex3 = os.path.isfile('./DTT/0'+filt+'/0'+days[x]+'_DAYS/NE/'+dt.strftime('%Y-%m-%d')+'.txt')
            ex4 = os.path.isfile('./DTT/0'+filt+'/0'+days[x]+'_DAYS/EE/'+dt.strftime('%Y-%m-%d')+'.txt')

            if ex1 and ex2 and ex3 and ex4:
                file = open('./DTT/0'+filt+'/0'+days[x]+'_DAYS/NN/'+dt.strftime('%Y-%m-%d')+'.txt',"r")
                line1 = file.readlines(); file.close();
                file = open('./DTT/0'+filt+'/0'+days[x]+'_DAYS/EN/'+dt.strftime('%Y-%m-%d')+'.txt',"r")
                line2 = file.readlines(); file.close();
                file = open('./DTT/0'+filt+'/0'+days[x]+'_DAYS/NE/'+dt.strftime('%Y-%m-%d')+'.txt',"r")
                line3 = file.readlines(); file.close();
                file = open('./DTT/0'+filt+'/0'+days[x]+'_DAYS/EE/'+dt.strftime('%Y-%m-%d')+'.txt',"r")
                line4 = file.readlines(); file.close();

                for p in range(0,len(sta_pairs)):
                    pair = sta_pairs[p]; az = azs[p];
                    p1 = False; p2 = False; p3 = False; p4 = False;

                    for k in range(1,len(line1)):
                        Date,Pairs,M,EM,A,EA,M0,EM0 = line1[k].split(",")
                        if Pairs.strip() == pair:
                            NNR = float(M)*math.cos(az)*math.cos(az)
                            NNT = float(M)*math.sin(az)*math.sin(az)
                            p1 = True;

                    for k in range(1,len(line2)):
                        Date,Pairs,M,EM,A,EA,M0,EM0 = line2[k].split(",")
                        if Pairs.strip() == pair:
                            EN = float(M)*math.cos(az)*math.sin(az)
                            p2 = True;

                    for k in range(1,len(line3)):
                        Date,Pairs,M,EM,A,EA,M0,EM0 = line3[k].split(",")
                        if Pairs.strip() == pair:
                            NE = float(M)*math.cos(az)*math.sin(az)
                            p3 = True;

                    for k in range(1,len(line4)):
                        Date,Pairs,M,EM,A,EA,M0,EM0 = line4[k].split(",")
                        if Pairs.strip() == pair:
                            EER = float(M)*math.sin(az)*math.sin(az)
                            EET = float(M)*math.cos(az)*math.cos(az)
                            p4 = True;

                    if p1 and p2 and p3 and p4:
                        assistR.append(NNR-EN-NE+EER)
                        assistT.append(NNT+EN+NE+EET)

            RR_pts.append(len(assistR)); RR_sigma[j] = -100*(mad(assistR)); RR_med[j] = -100*np.median(assistR)
            TT_pts.append(len(assistT)); TT_sigma[j] = -100*(mad(assistT)); TT_med[j] = -100*np.median(assistT)
            j += 1

        RR_sigma, RR_med, RR_pts = map(np.asarray, [RR_sigma, RR_med, RR_pts])
        TT_sigma, TT_med, TT_pts = map(np.asarray, [TT_sigma, TT_med, TT_pts])
        
        # Write data
        f1 = open('./results/f'+filt+'_m'+days[x]+'_all.txt', 'w')
        f1.write('time;ZZ_med;ZZ_sigma;ZZ_pts;RR_med;RR_sigma;RR_pts;TT_med;TT_sigma;TT_pts\n')

        for k in range(0,len(ZZ_med)):
            f1.write(str(ZZ_t[k])+';')
            f1.write(str(ZZ_med[k])+';'+str(ZZ_sigma[k])+';'+str(ZZ_pts[k])+';')
            f1.write(str(RR_med[k])+';'+str(RR_sigma[k])+';'+str(RR_pts[k])+';')
            f1.write(str(TT_med[k])+';'+str(TT_sigma[k])+';'+str(TT_pts[k])+'\n')

        f1.close()