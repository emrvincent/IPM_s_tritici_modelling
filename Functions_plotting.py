import matplotlib.pyplot as plt
from Functions_base import Temerge, T87, t_growing, T61, Amax
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from Functions_analysis import Y
import numpy as np
from matplotlib.ticker import PercentFormatter
from Functions_base import ic_base
from Functions_IPM import run_year_n,run_year_f

# Colours - green-beige
c1 = '#DACAA0'
c2 = '#B59441'
c3 = '#887D3D'
c4 = '#5A6638'
c5 = '#2D331C'

# Function to highlight grain forming period on input axis
def plot_grain_forming(ax,text=True):
    
    ax.axvspan(T61, T87, alpha=0.3, color=c1)
    if text == True:
        ax.text(2100,ax.get_ylim()[1]*0.92,"Grain forming",fontsize='9',c="#524C3C",weight='bold')
    return ax

# ###########################################
# ## Yield bar plot
# ###########################################
# def yield_barplot(ax,Ys,labels,title):
    
#     if len(Ys) != len(labels):
#         raise Exception("Not the right number of labels")
        
#     cols = [c3,c4,c5,c5,c5]
        
#     for i in range(len(Ys)):
#         ax.bar(range(i,i+1),Ys[i],color=cols[i],alpha=0.8)
        
#         tc = 'k'
#         ax.text(i, Ys[i] - 0.07,"%.3f" % round(Ys[i],3),fontsize='9',horizontalalignment='center',c=tc)
        
#     ax.set_title("Relative yield - "+title,fontsize='11')
#     ax.set_xticks(ticks = range(len(Ys)),labels = labels)
#     ax.set_ylim([0,1])
#     return ax

###########################################
## Infection prevalence figures for individual interventions
########################################### 
def plot_one_intervention(pop,labels,title,severity=1,fill=True):
    fig,ax1 = plt.subplots(1,1,figsize = (4.5,3))
    
    if len(pop) not in [3,4]:
        raise Exception("Need either 3 or 4 populations")
    if len(pop) != len(labels):
        raise Exception("Not the right number of labels")

    ls = ['--','-',':']
    cols = [c4,c5,c4]
        
    # Plot comparator no control and fungicide
    popf = run_year_f(ic_base,severity=severity)
    popn = run_year_n(ic_base,severity=severity)
    ax1.plot(t_growing,100*popf[:,2]/np.sum(popf[:,:5],axis=1),c=c2,label = "Fungicide",linewidth=1.5,marker='o',markevery=100,markersize=7)
    ax1.plot(t_growing,100*popn[:,2]/np.sum(popn[:,:5],axis=1),c=c2,label = "No control",linewidth=1.5,marker='^',markevery=100,markersize=7)
    
    # Plot the infection curves
    for i in range(len(pop)):
        tot = np.sum(pop[i][:,:5],axis=1)
        ax1.plot(t_growing,100*pop[i][:,2]/tot,c=cols[i],label = labels[i],linewidth=2.5,linestyle=ls[i])
    
    # Fill, if it's one with variable outcomes
    if fill == True:
        ax1.fill_between(t_growing, 100*pop[0][:,2]/np.sum(pop[0][:,:5],axis=1), 100*pop[-1][:,2]/np.sum(pop[-1][:,:5],axis=1), color=c4,alpha = 0.4,label="All possible\noutcomes")
    elif fill == False:
        None
    else:
        raise Exception("Fill type needs to be true or false")
    
    # Aesthetics
    newhandles,newlabels = ax1.get_legend_handles_labels()
    
    if fill == True:
        newhandles = newhandles[2:-1]
        newlabels = newlabels[2:-1]
    else:
        newhandles = newhandles[2:]
        newlabels = newlabels[2:]
#     newhandles = newhandles[2:]+newhandles[:2]
#     newlabels = newlabels[2:]+newlabels[:2]

    ax1.legend(handles = newhandles,labels=newlabels,loc= "upper left",fontsize = '9',handlelength=2.9)
    
    xmin = Temerge
    xmax = T87
    ax1.set_xlim([xmin, xmax])
    if (np.max(100*popn[:,2]/np.sum(popn[:,:5],axis=1)) < 20):
        ax1.set_ylim([0,25])
    else:
        ax1.set_ylim([0,30])
        
    # Plot grain forming
    ax1 = plot_grain_forming(ax1)
            
    ax1.set_title("Infection prevalence - " + title,fontsize='11')
    ax1.set_xlabel("Time (degree-days)")
    
    ax1.yaxis.set_major_formatter(PercentFormatter(decimals=0))
    
    # Print the yields
    Ys = [Y(i) for i in pop]
    print("Yields" + "\n" + str(Ys))
    
    # Print the peak infection prevalence and time
    peaks = [0]*len(pop)
    peaktimes = [0]*len(pop)
    for i in range(len(pop)):
        tot = np.sum(pop[i][:,:5],axis=1)
        I_perc = pop[i][:,2]/tot
        peaks[i] = np.max(I_perc)
        peaktimes[i] = np.argmax(I_perc) + Temerge
    print("Peak percent infection" + "\n" + str(peaks))
#     print("Peak infection time" + "\n" + str(peaktimes))
    
    return fig,ax1

# ###########################################
# ## Infection prevalence and yield for the single-field scenario
# ########################################### 
# def plot_one_field(pop,labels,title):
#     fig1,(ax1,ax12) = plt.subplots(1,2,figsize = (7,3))
    
#     if len(pop) not in [5]:
#         raise Exception("Need 3 populations, +2 variations for IPM")
#     if len(pop) != len(labels):
#         raise Exception("Not the right number of labels")
        
#     ls = ['-','--',':',':',':']
#     cols = [c3,c4,c5,c5,c5]
    
#     # Plot the infection curves
#     for i in range(len(pop)):
#         tot = np.sum(pop[i][:,:5],axis=1)
#         ax1.plot(t_growing,100*pop[i][:,2]/tot,c=cols[i],label = labels[i],linewidth=2.5,linestyle=ls[i])

#     ax1.legend(loc= "upper left",fontsize = '9')
    
#     xmin = Temerge
#     xmax = T87
#     ax1.set_xlim([xmin, xmax])
#     ax1.set_ylim([0,25])
#     ax1.set_title("Infection prevalence - " + title,fontsize='11')
#     ax1.set_xlabel("Time (degree-days)")
#     ax1.yaxis.set_major_formatter(PercentFormatter(decimals=0))
    
#     ax1 = plot_grain_forming(ax1)
    
#     # Plot the healthy curves
#     for i in range(len(pop)):
#         ax12.plot(t_growing,pop[i][:,0],c=cols[i],label = labels[i],linewidth=2.5,linestyle=ls[i])
  
#     xmin = Temerge
#     xmax = T87
#     ax12.set_xlim([xmin, xmax])
#     ax12.set_ylim([0,Amax])
#     ax12.set_title("Area healthy leaf tissue - " + title,fontsize='11')
#     ax12.set_xlabel("Time (degree-days)")
    
#     ax12 = plot_grain_forming(ax12)

    
#     # Plot the yields
#     fig2,ax2 = plt.subplots(1,1,figsize = (7,3))
#     Ys = [Y(i) for i in pop]
#     ax2 = yield_barplot(ax2,Ys,labels,title)
    
#     # Print the peak infection prevalence and time
#     peaks = [0]*len(pop)
#     peaktimes = [0]*len(pop)
#     for i in range(len(pop)):
#         tot = np.sum(pop[i][:,:5],axis=1)
#         I_perc = pop[i][:,2]/tot
#         peaks[i] = np.max(I_perc)
#         peaktimes[i] = np.argmax(I_perc) + Temerge
#     print("Peak percent infection" + "\n" + str(peaks))
#     print("Peak infection time" + "\n" + str(peaktimes))
    
#     return fig1,ax1,fig2,ax2

###########################################
## Plot SEI for one farm (Supp Info figure)
########################################### 
def plot_one_farm_type(pop, strategy_code):
    if strategy_code == "f":
        strategy = "Fungicide"
    elif strategy_code == "i":
        strategy = "IPM"
    elif strategy_code == "n":
        strategy = "No control"
    else:
        raise Exception("Input \"f\" or \"i\"")
        
    fig,ax = plt.subplots(1,figsize = (5,3))
    xmin = Temerge
    xmax = T87
    S = pop[:,0]
    E = pop[:,1]
    I = pop[:,2]
    ax.plot(t_growing,S,c=c5,label = "S",linewidth=2,linestyle=':')
    ax.plot(t_growing,E,c=c4,label = "E",linewidth=2,linestyle='-')
    ax.plot(t_growing,I,c=c2,label = "I",linewidth=2,linestyle='--')
    
#     # Plot inset
#     if strategy_code in ['f','i']:
#         inset_ax = inset_axes(ax,width='24%',height=0.6)
#         inset_ax.plot(t_growing,S,c='tab:green',label = "S",linewidth=2,linestyle=':')
#         inset_ax.plot(t_growing,E,c='tab:orange',label = "E",linewidth=2,linestyle='-')
#         inset_ax.plot(t_growing,I,c='tab:red',label = "I",linewidth=2,linestyle='--')
#         inset_ax.set_xlim([2500, xmax])
#         inset_ax.set_ylim([0,0.15])
#         inset_ax = plot_grain_forming(inset_ax,False)
    
    if strategy_code == 'n':
        ax.legend(loc = "upper right",fontsize = '9')
    
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([0,4.5])
    ax.set_title(strategy + " regime")
    ax.set_ylabel("Leaf area")
    ax.set_xlabel("Time (degree-days)")
    
    ax = plot_grain_forming(ax)
    
    return fig,ax