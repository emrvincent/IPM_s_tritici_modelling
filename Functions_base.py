#!/usr/bin/env python
# coding: utf-8

import math
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import minimize

###########################################
## Variables
###########################################

# Plant growth stages
Temerge = 1396
T31 = Temerge
T32 = 1495
T39 = 1653
T61 = 1891
T87 = 2567

# Variety mixture times
# T71 is 20th June in WGG
T71 = 2000
# We linearly interpolate between the closest ones in WGG
T69 = int(T61 + np.round((69-61)*(T71-T61)/(71-61)))
T75 = int(T71 + np.round((75-71)*(T87-T71)/(87-71)))

# Biocontrol times
# GS63 is 48 hours after GS61 - I found this one manually from the excel document
T63 = 1914

# Plant growth parameters
G = 8.2e-3
Amax = 4.438

# Disease system parameters
beta_base = 2.11e-2
gamma = 1/350
mu = 1/600
v = 8.97e-3

# Initial conditions
S0 = 0.05
psi = 1.44e-2
ic_base = [S0, 0, 0, 0, 0, psi]
ic_twofield = [S0, 0, 0, 0, 0, psi, S0, 0, 0, 0, 0, psi]

# Growing period
t_growing = np.arange(Temerge,T87,1)

# Fungicide parameters
omega_base = 1
theta_base = 7.677e+00
delta_base = 1.860e-02

###########################################
## Functions
###########################################

# Senescence function
def sigma(t):
    tau_s = 2.8e-3
    psi_s = 0.704
    omega_s = 0.314
    
    if t < T61:
        return 0
    else:
        return tau_s*(t-T61)/(T87-T61) + psi_s*math.exp(-omega_s*(T87-t))

# Plant growth function (basic)
def g_1D(A,t):
    if A == 0:
        return 0
    if t < Temerge:
        return 0
    else:
        return G*(Amax-A)
    
# Plant growth function (allows the proportion of the field type to be varied)
def g_2D(A,t,prop):
    if A == 0:
        return 0
    if t < Temerge:
        return 0
    else:
        return G*(prop*Amax-A)