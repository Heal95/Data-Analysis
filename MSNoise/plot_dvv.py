#!/usr/bin/env python
# coding: utf-8

# Plots dv/v data

# Import modules
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import os
from datetime import timedelta, date
import pylab
from obspy.geodetics.base import gps2dist_azimuth
import math

# Read dv/v files
def read_dvv(fn):
    with open(fn, 'r') as f1:
        lines = f1.readlines()
    time = []
    ZZ_sigma, ZZ_med, ZZ_pts = [], [], []
    RR_sigma, RR_med, RR_pts = [], [], []
    TT_sigma, TT_med, TT_pts = [], [], []
    for line in lines[1:]:
        t, z1, z2, z3, r1, r2, r3, t1, t2, t3 = line.split(';')
        time.append(t)
        ZZ_med.append(float(z1)); ZZ_sigma.append(float(z2)); ZZ_pts.append(float(z3))
        RR_med.append(float(r1)); RR_sigma.append(float(r2)); RR_pts.append(float(r3))
        TT_med.append(float(t1)); TT_sigma.append(float(t2)); TT_pts.append(float(t3))

    ZZ_sigma, ZZ_med, ZZ_pts = map(np.asarray, [ZZ_sigma, ZZ_med, ZZ_pts])
    RR_sigma, RR_med, RR_pts = map(np.asarray, [RR_sigma, RR_med, RR_pts])
    TT_sigma, TT_med, TT_pts = map(np.asarray, [TT_sigma, TT_med, TT_pts])

    sigma = [ZZ_sigma, RR_sigma, TT_sigma]
    median = [ZZ_med, RR_med, TT_med]
    points = [ZZ_pts, RR_pts, TT_pts]
    
    return time, median, sigma, points

# List of data to be considered
# Stacks
days = ['30', '60', '90']
# Filters
filters = ['0']
# Axis labels (need to match input files)
labels = ['ZZ','RR','TT']

# Plot data for each component, stack and filter
for filt in filters:
    for d in days:
        fn = './results/f'+filt+'_m'+d+'_all.txt'
        time, median, sigma, points = read_dvv(fn)

        fig = plt.figure(figsize=(20, 15))
        plt.rcParams.update({'font.size': 18})
        
        for k in range(1,4):
            plt.subplot(3,1,k)
            if k==1:
                plt.title(d+' day stack, filter f'+filt)
            plt.grid(True)
            plt.scatter(time, median[k-1], s=20, c=points[k-1], cmap='gist_rainbow_r',label=labels[k-1])
            plt.fill_between(time, 1.0*(median[k-1]-sigma[k-1]), 1.0*(median[k-1]+sigma[k-1]), zorder=-1, alpha=0.6, color ='gray', lw=0)
            axes = plt.gca()
            axes.set_ylim([-0.2,0.2])
            start, end = axes.get_xlim()
            axes.xaxis.set_ticks(np.arange(start-183, end+10, 2*365.25))
            axes.set_xlim([-200,None])
            plt.legend(loc='upper right')
            cbar = plt.colorbar()
            cbar.ax.set_ylabel('Station pairs', rotation=270, labelpad=30)
            plt.ylabel('dv/v [%]')

        plt.tight_layout()
        plt.show()
        fig.savefig('./plots/f'+filt+'_m'+d+'_all.png')