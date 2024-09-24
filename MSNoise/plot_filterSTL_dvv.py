#!/usr/bin/env python
# coding: utf-8

# Plots low-pass filtered dv/v data and its STL decomposition

# Import modules
import os
import math
import obspy
import statsmodels
import numpy as np
import pandas as pd
import obspy.signal
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from statsmodels.tsa.seasonal import STL
from obspy import Stream, Trace, UTCDateTime

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

# Convert dataset to Obspy trace object
def data_to_trace(data,delta,net,sta,label,starttime):
    trace = Trace(data=data)
    trace.stats.delta = delta
    trace.stats.network = net
    trace.stats.station = sta
    trace.stats.location = ""
    trace.stats.channel = label
    trace.stats.starttime = starttime
    tr = trace
    
    return tr

# List of data to be considered
# Stacks
stacks = [30, 60, 90]
#Filters
filters = [0]
# Axis labels (need to match input files)
labels = ['ZZ','RR','TT']
# Period for low-pass filtering
T = 365.25
freqmin = 1./(1.0*3600*24*0.5*T)

# Plot data for each component, stack and filter
for d in stacks:
    for filt in filters:
        fn = './results/f'+str(filt)+'_m'+str(d)+'_all.txt'
        time, median, sigma, points = read_dvv(fn)
        starttime = UTCDateTime(2009, 1, 1, 0, 0, 0) + timedelta(days=d-1)
        
        # Plot original and low-pass filtered data
        fig = plt.figure(figsize=(20,30))
        plt.rcParams.update({'font.size': 18})
        
        for k in range(0,3):
            comp_med = median[k]; 
            comp_med = comp_med[(d-1):]

            for j in range(0,len(comp_med)):
                if str(comp_med[j])=='nan':
                    comp_med[j] = (comp_med[j-1]+comp_med[j+1])/2

            tr = data_to_trace(comp_med, 1.0*3600*24, "MSN", "ALL", labels[k], starttime)
            tr_comp = tr.detrend("simple")
        
            tr_compf = tr_comp.copy()
            tr_compf.filter("lowpass", freq = freqmin, corners=4, zerophase=True)

            ax = fig.add_subplot(3, 1, k+1)
            if k==0:
                plt.title(str(d)+' day stack, filter f'+str(filt))
            ax.plot(tr_comp.times("matplotlib"), tr_comp.data, "k",lw=0.5, label=labels[k])
            ax.plot(tr_compf.times("matplotlib"), tr_compf.data, "tab:orange",lw=4, label='LP filtered '+labels[k])
            plt.ylabel('dv/v [%]')
            ax.xaxis_date()
            fig.autofmt_xdate()
            plt.legend(loc='upper right')
            plt.grid()

        plt.show()
        fig.savefig('./plots/f'+str(filt)+'_m'+str(d)+'.jpg')
        plt.close()

        # Plot decomposition of dataset
        fig=plt.figure(figsize=(60, 20))
        plt.rcParams.update({'font.size': 12})

        for k in range(0,3):
            comp_med = median[k]; 
            comp_med = comp_med[(d-1):]

            for j in range(0,len(comp_med)):
                if str(comp_med[j])=='nan':
                    comp_med[j] = (comp_med[j-1]+comp_med[j+1])/2

            tr = data_to_trace(comp_med, 1.0*3600*24, "MSN", "ALL", labels[k], starttime)
            tr_comp = tr.detrend("simple")
        
            # Convert the signal to a pandas dataframe
            df_comp = pd.DataFrame(tr_comp.data, index=pd.date_range(tr_comp.stats.starttime.strftime("%Y-%m-%d"),periods=len(tr_comp.data), freq='D'), columns=['signal'])
       
            # Decompose the signal using the STL method
            result_comp = STL(df_comp,period=365,robust=True).fit()

            # Plot the original signal and the extracted seasonal, trend, and residual components
            ax1 = fig.add_subplot(3, 4, 4*k+1)
            if k==0:
                plt.title('Original signal')
            ax1.plot(df_comp, label=labels[k], c="k", lw=0.5)
            plt.ylabel('dv/v [%]')
            ax1.xaxis_date()
            fig.autofmt_xdate()
            plt.legend(loc='upper right')

            ax2 = fig.add_subplot(3, 4, 4*k+2)
            if k==0:
                plt.title('Seasonal component')
            ax2.plot(result_comp.seasonal, c="b", lw=0.5)
            ax2.xaxis_date()
            fig.autofmt_xdate()

            ax3 = fig.add_subplot(3, 4, 4*k+3)
            if k==0:
                plt.title('Trend component')
            ax3.plot(result_comp.trend, c="r", lw=2.5)
            ax3.xaxis_date()
            fig.autofmt_xdate()

            ax4 = fig.add_subplot(3, 4, 4*k+4)
            if k==0:
                plt.title('Residual component')
            ax4.plot(result_comp.resid, c="tab:purple", lw=0.5)
            ax4.xaxis_date()
            fig.autofmt_xdate()

        plt.show()
        fig.savefig('./plots/f'+str(filt)+'_m'+str(d)+'_STL.jpg')
        plt.close()