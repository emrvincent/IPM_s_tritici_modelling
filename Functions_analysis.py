from Functions_base import T61,T87,Temerge
import numpy as np

nodis_Y = np.load("baseY.npy")

###########################################
## Compute yield
###########################################
# Yield computation function
def integral(S,E):
    # Define the grain forming period
    grain_period = np.arange(T61,T87,1)

    # Find indices of grain period
    t1 = T61-Temerge
    t2 = T87-Temerge

    # Compute yield
    temp_Y = np.trapz((S+E)[t1:t2], grain_period)
    return temp_Y
    
# Relative yield
def Y(disease_pop):
    # Validity check
    if disease_pop.shape[1] != 6:
        raise Exception("This needs to be size 6")
    
    # Yield
    S,E = disease_pop[:,:2].T
    this_Y = integral(S, E)

    return this_Y/nodis_Y

# ###########################################
# ## Other potential analysis functions
# ###########################################

# # Integral under infection
# def integral_I(disease_pop):
#     if disease_pop.shape[1] != 6:
#         raise Exception("This needs to be size 6")
#     return np.trapz(disease_pop[:,2])

# # % of dead crop
# def percent_D(disease_pop):
#     if disease_pop.shape[1] != 6:
#         raise Exception("This needs to be size 6")
#     return disease_pop[:,-2]/np.sum(disease_pop[:,:-1],axis=1)

# # Peak time
# def peak_time(disease_pop):
#     if disease_pop.shape[1] != 6:
#         raise Exception("This needs to be size 6")
#     return np.argmax(disease_pop[:,2]) + Temerge
