#!/usr/bin/env python

#module Loc_eq_f
#coding=utf-8

# Importing modules
import numpy as np
import math
import os

# Model info (Balkan earth model)
h1 = 30; vp1 = 5.40; vs1 = 3.05
h2 = 16; vp2 = 6.66; vs2 = 3.83
vp3 = 8.00; vs3 = 4.60

# Functions
def delta(phi, phi_s, lamb, lamb_s):
    x = np.sin(math.radians(phi))*np.sin(math.radians(phi_s)) + np.cos(math.radians(phi))*np.cos(math.radians(phi_s))*np.cos(math.radians(lamb-lamb_s))
    y = math.degrees(math.acos(x))
    return y

def azimuth(phi_s, phi, lamb, lamb_s, delta):
    x = np.sin(math.radians(lamb_s-lamb))*np.cos(math.radians(phi_s))/np.sin(math.radians(delta))
    if (phi_s-phi)<0:
        y = 180-math.degrees(math.asin(x))
    elif (phi_s-phi)>0 and (lamb_s-lamb)<0:
        y = 360+math.degrees(math.asin(x))
    elif (phi_s-phi)>0 and (lamb_s-lamb)>0:
        y = math.degrees(math.asin(x))
    return y

def tg(v1, delta, D): 
    y = (1/v1) * math.sqrt((111.2*delta)**2 + D**2)
    return y

def tn(v1, v2, v3, h1, h2, delta, D):
    x1 = (111.2*delta)/v3
    x2 = (2*h1 - D)*(math.sqrt(v3**2-v1**2)/(v1*v3))
    x3 = 2*h2*(math.sqrt(v3**2-v2**2)/(v2*v3))
    y = x1 + x2 + x3
    return y

def aPg(v1, delta, D, azimuth, phi): 
    x = ((delta*111.2)/v1) * (1/math.sqrt((111.2*delta)**2 + D**2))
    y = -x * np.sin(math.radians(azimuth)) * np.cos(math.radians(phi))
    return y

def aPn(v3, azimuth, phi):
    y = - np.sin(math.radians(azimuth))*np.cos(math.radians(phi))/v3
    return y

def bPg(v1, delta, D, azimuth):
    x = ((delta*111.2)/v1) * (1/math.sqrt((111.2*delta)**2 + D**2))
    y = -x * np.cos(math.radians(azimuth))
    return y

def bPn(v3, azimuth):
    y = - np.cos(math.radians(azimuth))/v3
    return y

def cPg(v1, delta, D):
    y = (D/v1) * (1/math.sqrt((111.2*delta)**2 + D**2))
    return y

def cPn(v2, v3):
    y = - math.sqrt(v3**2-v2**2)/(v2*v3)
    return y