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
Temerge = 1550
T31 = Temerge
T32 = 1668
T39 = 1912
T61 = 2255
T87 = 3104

# Plant growth parameters
R = 1.26e-2
k = 4.1

# Disease system parameters
beta_base = 2.08e-2
gamma = 1/266
mu = 1/456
v = 8.5e-3

# Initial conditions
S0 = 0.05
psi = 1.09e-2
ic_base = [S0, 0, 0, 0, 0, psi]
ic_twofield = [S0, 0, 0, 0, 0, psi, S0, 0, 0, 0, 0, psi]

# Growing period
t_growing = np.arange(Temerge,T87,1)

# Fungicide parameters
omega_base = 1
theta_base = 9.6
delta_base = 1.11e-2

###########################################
## Functions
###########################################

# Senescence function
def sigma(t):
    if t < T61:
        return 0
    else:
        return 0.005*(t-T61)/(T87-T61) + 0.1*math.exp(-0.02*(T87-t))

# Plant growth function (basic)
def g_1D(A,t):
    if A == 0:
        return 0
    if t < Temerge:
        return 0
    else:
        return R*(k-A)
    
# Plant growth function (allows the proportion of the field type to be varied)
def g_2D(A,t,prop):
    if A == 0:
        return 0
    if t < Temerge:
        return 0
    else:
        return R*(prop*k-A)