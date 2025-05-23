import numpy as np
from scipy.integrate import odeint
import math

# For ODEs
from Functions_base import beta_base, Temerge, gamma, g_1D, sigma, mu, v, T31, T39, T32, g_2D, T87, omega_base, theta_base, delta_base, t_growing

# Load ascospore data
prop_asc = np.load("prop_asc.npy")
growth_delay = np.load("dd_delay.npy")

###########################################
## No control system
###########################################
# No control ODEs
def dPop_n(x, t, severity):
    beta = beta_base*severity
    
    S,E,I,R,D,P = x
    A = S + E + I + R + D

    # Uncontrolled fields
    if t > Temerge:
        
        # Set disease development rate (force of latency)
        FOL = gamma*E
        
        # Set transmission rate (force of infection)
        FOI = beta*(I + P)/A
        
        dS = g_1D(A,t) - sigma(t)*S - S*FOI
        dE = S*FOI - sigma(t)*E - FOL
        dI = FOL - mu*I
        dR = sigma(t)*(S+E)
        dD = mu*I
        dP = -v*P
        
    else:
        dS = 0
        dE = 0
        dI = 0
        dR = 0
        dD = 0
        dP = 0

    return [dS,dE,dI,dR,dD,dP]

# Run no control for a year
def run_year_n(ic,severity = 1):
    pop = odeint(func = dPop_n, y0 = ic, t = t_growing, args = (severity,))
    return pop
    
###########################################
## Fungicide system
###########################################
# Effect of fungicide
def eps(t,omega,theta,delta):
    if t < T32:
        C = 0
    
    elif t < T39:
        C = math.exp(-delta*(t-T32))
    
    else:
        C = math.exp(-delta*(t-T32)) + math.exp(-delta*(t-T39)) 
    
    E = omega*(1-math.exp(-theta*C))

    return E

# Fungicide ODEs
def dPop_f(x, t,omega,theta,delta,severity):

    beta = beta_base*severity
    
    S,E,I,R,D,P = x
    A = S + E + I + R + D

    # Fungicide fields
    if t > Temerge:
        
        # Set disease development rate (force of latency)
        FOL = gamma*E*(1 - eps(t,omega,theta,delta))
        
        # Set transmission rate (force of infection)
        FOI = beta*(I + P)*(1 - eps(t,omega,theta,delta))/A
        
        dS = g_1D(A,t) - sigma(t)*S - S*FOI
        dE = S*FOI - sigma(t)*E - FOL
        dI = FOL - mu*I
        dR = sigma(t)*(S+E)
        dD = mu*I
        dP = -v*P
        
    else:
        dS = 0
        dE = 0
        dI = 0
        dR = 0
        dD = 0
        dP = 0

    return [dS,dE,dI,dR,dD,dP]

# Run fungicide for a year
def run_year_f(ic, omega = omega_base, theta = theta_base, delta = delta_base, severity = 1):
    pop = odeint(func = dPop_f, y0 = ic, t = t_growing,args=(omega,theta,delta,severity))
    return pop

###########################################
## IPM system
###########################################
# IPM ODEs
def dPop_i(x,t,r_beta,days_late,biocontrol,severity):
    
    beta = beta_base*severity
    
    S,E,I,R,D,P = x
    A = S + E + I + R + D
    delay = growth_delay[days_late]
    
    # IPM fields
    # Set sowing date dependant time delay
    if t > Temerge + delay:
        
        # Set biocontrol dependant disease rate
        FOL = gamma*(1 - biocontrol(t,delay))
            
        # Set variety and biocontrol dependant transmission rate
        FOI = r_beta*beta*(I + P)/A

        dS = g_1D(A,t) - sigma(t)*S - r_beta*S*FOI*(1 - biocontrol(t,delay))
        dE = r_beta*S*FOI*(1 - biocontrol(t,delay)) - sigma(t)*E - FOL*E
        dI = FOL*E - mu*I
        dR = sigma(t)*(S+E)
        dD = mu*I
        dP = -v*P
        
    else:
        dS = 0
        dE = 0
        dI = 0
        dR = 0
        dD = 0
        dP = 0

    return [dS,dE,dI,dR,dD,dP]

# Run IPM for a year
def run_year_i(input_ic,r_beta,days_late,prop_debris,cleared_debris,biocontrol_dates,severity=1):
    ic = 1*input_ic
    this_biocontrol = biocontrol(biocontrol_dates)

    p = prop_debris
    q = cleared_debris

    # Apply debris clearing
    ic = ic_IPM(ic,p,q,days_late,six_permission = True)
    
    # Run simulation for 1 year
    pop = odeint(func = dPop_i, y0 = ic, t = t_growing, args = (r_beta,days_late,this_biocontrol,severity))
    return pop

# Effect of biocontrol
def biocontrol(biocontrol_dates):
    def biocontrol_sub(t,delay):
        # Load the biocontrol parameters fitted for in Fitting_4_biocontrol
        omega, theta, delta =  np.load("biocontrol_fit_results.npy")

        # Check that biocontrol dates are valid
        for i in biocontrol_dates:
            if i not in [31,39]:
                raise Exception("Can only have 31 or 39 as biocontrol dates")

        # Switch biocontrol on for the relevant dates
        if 31 in biocontrol_dates:
            d31 = 1
        else:
            d31 = 0

        if 39 in biocontrol_dates:
            d39 = 1
        else:
            d39 = 0

        # Fungicide concentration
        if t < T31 + delay:
            C = 0
        elif t < T39 + delay:
            C = d31*math.exp(-delta*(t-(T31+delay)))
        else:
            C = d31*math.exp(-delta*(t-(T31+delay))) + d39*math.exp(-delta*(t-(T39+delay)))

        return omega*(1-math.exp(-theta*C))
    return biocontrol_sub

# Apply IPM to initial condition for ascospores
def ic_IPM(input_ic,p,q,days_late,six_permission = False):
    ic = 1*input_ic
    if len(ic) != 12 and six_permission == False:
        raise Exception("Need the input to just be the 12 IPM entries")

    # Apply debris clearing
    ic[5] = p*q*ic[5] + (1-p)*ic[5]
    
    # Apply the loss in ascospores due to late planting
    ic[5] = ic[5]*prop_asc[days_late]
    return ic

###########################################
## IPM+fungicide spot treatment system
###########################################
# Effect of fungicide
def eps_spot(t,omega_spot,Tspot):
    theta = theta_base
    delta = delta_base
    
    if t < Tspot:
        C = 0
    
    else:
        C = math.exp(-delta*(t-Tspot))
    
    E = omega_spot*(1-math.exp(-theta*C))

    return E

# IPM+fungicide spot treatment ODEs
def dPop_fi(x,t,r_beta,days_late,biocontrol,omega_spot,Tspot,severity):
    
    beta = beta_base*severity
    
    S,E,I,R,D,P = x
    A = S + E + I + R + D
    delay = growth_delay[days_late]
    
    # IPM fields
    # Set sowing date dependant time delay
    if t > Temerge + delay:
        
        # Set biocontrol dependant disease rate
        FOL = gamma*(1 - biocontrol(t,delay))*(1-eps_spot(t,omega_spot,Tspot))
            
        # Set variety and biocontrol dependant transmission rate
        FOI = r_beta*beta*(I + P)*(1 - biocontrol(t,delay))*(1-eps_spot(t,omega_spot,Tspot))/A

        dS = g_1D(A,t) - sigma(t)*S - r_beta*S*FOI
        dE = r_beta*S*FOI - sigma(t)*E - FOL*E
        dI = FOL*E - mu*I
        dR = sigma(t)*(S+E)
        dD = mu*I
        dP = -v*P
        
    else:
        dS = 0
        dE = 0
        dI = 0
        dR = 0
        dD = 0
        dP = 0

    return [dS,dE,dI,dR,dD,dP]

# Run IPM+fungicide for a year
def run_year_fi(input_ic,r_beta,days_late,prop_debris,cleared_debris,biocontrol_dates,omega_spot,Tspot,severity=1):
    ic = 1*input_ic
    this_biocontrol = biocontrol(biocontrol_dates)

    p = prop_debris
    q = cleared_debris

    # Apply debris clearing
    ic = ic_IPM(ic,p,q,days_late,six_permission = True)
    
    # Run simulation for 1 year
    pop = odeint(func = dPop_fi, y0 = ic, t = t_growing, args = (r_beta,days_late,this_biocontrol,omega_spot,Tspot,severity))
    return pop

#####################################
## Cross-field ODE system
#####################################
# Disease system with IPM and fungicide field
def dPop_mix2(x,t,r_beta,delay,biocontrol,prop_i,omega,theta,delta,severity):
    # Note that the input initial conditions should already have been scaled
    
    beta = beta_base*severity
    
    # Set initial conditions
    S_F,E_F,I_F,R_F,D_F,P_F = x[6:]
    A_F = S_F + E_F + I_F + R_F + D_F
    
    S_I,E_I,I_I,R_I,D_I,P_I = x[:6]
    A_I = S_I + E_I + I_I + R_I + D_I
    
    # IPM fields
    # Set sowing date dependant time delay
    if t >= Temerge + delay:
        
        # Set biocontrol dependant disease rate
        FOL = gamma*(1 - biocontrol(t,delay))
            
        # Set variety and biocontrol dependant transmission rate
        FOI = beta*(r_beta*(I_I + P_I) + I_F + P_F)/(A_I+A_F)
        
        
        dS_I = g_2D(A_I,t,prop_i) - sigma(t)*S_I - r_beta*S_I*FOI*(1 - biocontrol(t,delay))
        dE_I = r_beta*S_I*FOI*(1 - biocontrol(t,delay)) - sigma(t)*E_I - FOL*E_I
        dI_I = FOL*E_I - mu*I_I
        dR_I = sigma(t)*(S_I+E_I)
        dD_I = mu*I_I
        dP_I = -v*P_I
        
    else:
        dS_I = 0
        dE_I = 0
        dI_I = 0
        dR_I = 0
        dD_I = 0
        dP_I = 0

    # Fungicide fields
    if t >= Temerge:
        
        # Set fungicide dependant disease rate
        FOL = gamma*(1 - eps(t,omega,theta,delta))
        
        # Set fungicide dependant FOI
        FOI = beta*(r_beta*(I_I + P_I) + I_F + P_F)/(A_I+A_F)
        
        dS_F = g_2D(A_F,t,1-prop_i) - sigma(t)*S_F - S_F*FOI*(1 - eps(t,omega,theta,delta))
        dE_F = S_F*FOI*(1 - eps(t,omega,theta,delta)) - sigma(t)*E_F - FOL*E_F
        dI_F = FOL*E_F - mu*I_F
        dR_F = sigma(t)*(S_F+E_F)
        dD_F = mu*I_F
        dP_F = -v*P_F
        
    else:
        dS_F = 0
        dE_F = 0
        dI_F = 0
        dR_F = 0
        dD_F = 0
        dP_F = 0

    return [dS_I,dE_I,dI_I,dR_I,dD_I,dP_I, dS_F,dE_F,dI_F,dR_F,dD_F,dP_F]

# Apply the proportion of IPM to the initial conditions
def ic_proportion(input_ic,prop_i):
    ic = 1*input_ic
    if len(ic) != 12:
        raise Exception("This function needs length 12 vector")
        
    if any(ic[1:5]) != 0 or any(ic[7:11])!= 0:
        raise Exception("Needs to be zero outside of S and P")
                        
    # Apply the relevant proportion of farms in each control method
    ic[0] = prop_i*ic[0]
    ic[5] = prop_i*ic[5]
    
    ic[6] = (1-prop_i)*ic[6]
    ic[11] = (1-prop_i)*ic[11]
    
    return ic

# Run cross-field scenario for 1 year
def run_year_mix2(input_ic,r_beta,days_late,prop_debris,cleared_debris,biocontrol_dates,prop_i,omega=omega_base,theta=theta_base,delta=delta_base,severity=1):
    
    if prop_i in [0,1]:
        raise Exception("Don't put prop_i=0,1 in here")
    ic = 1*np.array(input_ic)
    
    # Apply IPM: biocontrol and initial conditions
    this_biocontrol = biocontrol(biocontrol_dates)
    ic = ic_IPM(ic,prop_debris,cleared_debris,days_late)
    
    if prop_i == 0:
        ic[:6] = np.zeros(6)
    elif prop_i == 1:
        ic[6:] = np.zeros(6)
    
    # Proportion the whole system according to the IPM/fungicide ratio
    ic = ic_proportion(ic,prop_i)
    
    # Run first period before IPM fields are sown
    delay = growth_delay[days_late]
    ic0 = 1*ic
    ic0[:6] = np.zeros(6)
    t_0 = np.arange(Temerge,Temerge+delay+1,1)
    pop0 = odeint(func = dPop_mix2, y0 = ic0, t = t_0, args = (r_beta,delay,this_biocontrol,prop_i,omega,theta,delta,severity))
    
    # Run second period after IPM fields are sown
    ic1 = 1*pop0[-1]
    ic1[:6] = 1*ic[:6]
    t_1 = np.arange(Temerge+delay,T87,1)
    pop1 = odeint(func = dPop_mix2, y0 = ic1, t = t_1, args = (r_beta,delay,this_biocontrol,prop_i,omega,theta,delta,severity))
    
    # Concatenate
    pop = np.concatenate((pop0[:-1,:],pop1))
    t = np.concatenate((t_0[:-1],t_1))
    return t,pop

#####################################
## Cross-field IPM+fungicide ODE system
#####################################
# Disease system with IPM and fungicide field
def dPop_mix_fi(x,t,r_beta,delay,biocontrol,prop_i,omega_spot,Tspot,severity):
    # Note that the input initial conditions should already have been scaled
    
    beta = beta_base*severity
    omega = omega_base
    theta = theta_base
    delta = delta_base
    
    # Set initial conditions
    S_F,E_F,I_F,R_F,D_F,P_F = x[6:]
    A_F = S_F + E_F + I_F + R_F + D_F
    
    S_I,E_I,I_I,R_I,D_I,P_I = x[:6]
    A_I = S_I + E_I + I_I + R_I + D_I
    
    # IPM fields
    # Set sowing date dependant time delay
    if t >= Temerge + delay:
        
        # Set biocontrol dependant disease rate
        FOL = gamma*(1 - biocontrol(t,delay))*(1-eps_spot(t,omega_spot,Tspot))
            
        # Set variety and biocontrol dependant transmission rate
        FOI = beta*(r_beta*(I_I + P_I) + I_F + P_F)*(1-eps_spot(t,omega_spot,Tspot))/(A_I+A_F)
        
        
        dS_I = g_2D(A_I,t,prop_i) - sigma(t)*S_I - r_beta*S_I*FOI*(1 - biocontrol(t,delay))
        dE_I = r_beta*S_I*FOI*(1 - biocontrol(t,delay)) - sigma(t)*E_I - FOL*E_I
        dI_I = FOL*E_I - mu*I_I
        dR_I = sigma(t)*(S_I+E_I)
        dD_I = mu*I_I
        dP_I = -v*P_I
        
    else:
        dS_I = 0
        dE_I = 0
        dI_I = 0
        dR_I = 0
        dD_I = 0
        dP_I = 0

    # Fungicide fields
    if t >= Temerge:
        
        # Set fungicide dependant disease rate
        FOL = gamma*(1 - eps(t,omega,theta,delta))
        
        # Set fungicide dependant FOI
        FOI = beta*(r_beta*(I_I + P_I) + I_F + P_F)/(A_I+A_F)
        
        dS_F = g_2D(A_F,t,1-prop_i) - sigma(t)*S_F - S_F*FOI*(1 - eps(t,omega,theta,delta))
        dE_F = S_F*FOI*(1 - eps(t,omega,theta,delta)) - sigma(t)*E_F - FOL*E_F
        dI_F = FOL*E_F - mu*I_F
        dR_F = sigma(t)*(S_F+E_F)
        dD_F = mu*I_F
        dP_F = -v*P_F
        
    else:
        dS_F = 0
        dE_F = 0
        dI_F = 0
        dR_F = 0
        dD_F = 0
        dP_F = 0

    return [dS_I,dE_I,dI_I,dR_I,dD_I,dP_I, dS_F,dE_F,dI_F,dR_F,dD_F,dP_F]

# Run cross-field scenario for 1 year
def run_year_mix_fi(input_ic,r_beta,days_late,prop_debris,cleared_debris,biocontrol_dates,prop_i,omega_spot,Tspot,severity):
    
    if prop_i in [0,1]:
        raise Exception("Don't put prop_i=0,1 in here")
    ic = 1*np.array(input_ic)
    
    # Apply IPM: biocontrol and initial conditions
    this_biocontrol = biocontrol(biocontrol_dates)
    ic = ic_IPM(ic,prop_debris,cleared_debris,days_late)
    
    if prop_i == 0:
        ic[:6] = np.zeros(6)
    elif prop_i == 1:
        ic[6:] = np.zeros(6)
    
    # Proportion the whole system according to the IPM/fungicide ratio
    ic = ic_proportion(ic,prop_i)
    
    # Run first period before IPM fields are sown
    delay = growth_delay[days_late]
    ic0 = 1*ic
    ic0[:6] = np.zeros(6)
    t_0 = np.arange(Temerge,Temerge+delay+1,1)
    pop0 = odeint(func = dPop_mix_fi, y0 = ic0, t = t_0, args = (r_beta,delay,this_biocontrol,prop_i,omega_spot,Tspot,severity))
    
    # Run second period after IPM fields are sown
    ic1 = 1*pop0[-1]
    ic1[:6] = 1*ic[:6]
    t_1 = np.arange(Temerge+delay,T87,1)
    pop1 = odeint(func = dPop_mix_fi, y0 = ic1, t = t_1, args = (r_beta,delay,this_biocontrol,prop_i,omega_spot,Tspot,severity))
    
    # Concatenate
    pop = np.concatenate((pop0[:-1,:],pop1))
    t = np.concatenate((t_0[:-1],t_1))
    return t,pop
