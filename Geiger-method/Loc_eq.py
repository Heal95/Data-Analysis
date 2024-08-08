#!/usr/bin/env python
"""
Created on Fri Oct 22 09:45:45 2021

@author: Helena

A script to relocate earthquakes using Geiger method. Input files are bulletins from SANDI.

The script uses the Iterate module stored in the same working directory.

"""

# Import modules
import Iterate as it
import os
import numpy as np

# Input folder where Bulletins are
in_folder = './BULLETINS'
subdir = os.listdir(in_folder)

for file in subdir:
    Depth = 10 # Assumed hypocenter depth
    
    eq_hyp_time, phi_eq_pp, lamb_eq_pp, T, sta_phase = it.read_bulletin(in_folder, file)
    delta_sta, azimuth_sta, stations = it.sta_azimuth_delta(sta_phase, phi_eq_pp, lamb_eq_pp)
    eq_hyp_time_pp = it.fix_origin_time(sta_phase, T, eq_hyp_time, delta_sta, Depth)

    sigma = 0.001
    X, lamb_eq_pp, phi_eq_pp, Depth, eq_hyp_time_pp, i = it.Geiger_first_it(sta_phase, T, eq_hyp_time_pp, phi_eq_pp, lamb_eq_pp, delta_sta, azimuth_sta, Depth)
    eq_hyp_time_pp, phi_eq_pp, lamb_eq_pp, Depth = it.Geiger_it(X, sigma, i, stations, phi_eq_pp, lamb_eq_pp, sta_phase, T, eq_hyp_time_pp, Depth)
    
    print("Event:", file)
    mins = ((eq_hyp_time_pp[0]/3600)-int((eq_hyp_time_pp[0]/3600)))*60
    print(mins)
    secs = (mins-int(mins))*60
    print("Origin time:", file[-4:-2], "h", file[-2:], "min", round(secs,4), "s")
    print("Latitude:", phi_eq_pp[0])
    print("Longitude:", lamb_eq_pp[0])
    print("Depth:", Depth[0], "km")
    print("------------------------------")