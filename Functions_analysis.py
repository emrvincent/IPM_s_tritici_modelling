from Functions_base import T61,T87,ic_base,Temerge
from Functions_IPM import run_year_n
import numpy as np

###########################################
## Compute yield
###########################################

# Yield
def Y(disease_pop):
    # Validity check
    if disease_pop.shape[1] != 6:
        raise Exception("This needs to be size 6")

    # Define the grain forming period
    grain_period = np.arange(T61,T87,1)
    
    # Find indices of grain period
    t1 = T61-Temerge
    t2 = T87-Temerge
    
    # Yield computation function
    def integral(S,E):
        # Compute yield
        temp_Y = np.trapz((S+E)[t1:t2], grain_period)
        return temp_Y
    
    # Yield
    S,E = disease_pop[:,:2].T
    this_Y = integral(S, E)
    
    # Uninfected yield for denominator
    ic = 1*ic_base
    ic[-1] = 0
    disease_free_pop = run_year_n(ic)
    S_no_dis = disease_free_pop[:,0]
    nodis_Y = integral(S_no_dis,0)
    
    return this_Y/nodis_Y

###########################################
## Other potential analysis functions
###########################################

# Integral under infection
def integral_I(disease_pop):
    if disease_pop.shape[1] != 6:
        raise Exception("This needs to be size 6")
    return np.trapz(disease_pop[:,2])

# % of dead crop
def percent_D(disease_pop):
    if disease_pop.shape[1] != 6:
        raise Exception("This needs to be size 6")
    return disease_pop[:,-2]/np.sum(disease_pop[:,:-1],axis=1)

# Peak time
def peak_time(disease_pop):
    if disease_pop.shape[1] != 6:
        raise Exception("This needs to be size 6")
    return np.argmax(disease_pop[:,2]) + Temerge