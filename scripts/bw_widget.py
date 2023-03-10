import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from leiva2020 import Two_layer
from dqdv_proc import DQDV
from scipy.signal import chirp, find_peaks, peak_widths
import os
import pandas as pd

R = 8.31 # Molar gas constant.
delta_eta = 0.02 # spacing on the axes
eta_1_init =0.7 # initial value of eta
F = 96.485 # Faraday constant, kC per mole

# Initialisation of arguments.
arg_dict = {}
used_args= ['E0','delE','g2','alpha4','beta4','T','M','L']
unused_args = ['g1','alpha3','beta3','alpha1','beta1']

entropy_path = '../entropy_data/'
figures_path = '../figures/'
galvan_path = '../galvanostatic/'
experimental_dir = entropy_path +'entropydischarge_hc_00/'
mass = 0.9*0.004203 # Mass of carbon active material, in g, 90% loading

for arg in unused_args:  # Pad with zeros
    arg_dict[arg] = [0]

arg_dict['E0'] = [-0.25]
arg_dict['delE'] = [-0.25]
arg_dict['g2'] = [-0.020]
arg_dict['alpha4'] = [-1.0]
arg_dict['beta4'] = [2.0]
arg_dict['T'] = [298.0]
arg_dict['M'] = 120
arg_dict['L'] = [0.250]
arg_dict['label'] = 'L'
initial_cap = 250
T = arg_dict['T'][0] # T as float, in K

two_layer = Two_layer(arg_dict) # Instantiate the class
df_bw = two_layer.solution() # Run all thermodynamics
suffix = str(format(arg_dict[arg_dict['label']][0],'.3f'))
print(list(df_bw))

def half_widths(dQdV):        
    peaks, _ = find_peaks(dQdV)
    results_half = peak_widths(self.dQdV[g], peaks, rel_height=0.5)
    left_index = int(results_half[2])
    right_index = int(results_half[3])
    E_g = self.E[g] # Change to read the simulated electrode potential
    half_width = -(E_g[right_index] - E_g[left_index]) * 1000
    self.peak_half_widths.append(half_width)
#            print('g='+ str(g)+ ', width ='+ str(results_half[0]) + ', width_heights='+ str(results_half[1]) + ', left_lps=' + str(results_half[2]) + ', right_lps=' + str(results_half[3]))
    print('width ='+ str(half_width))

def Sconfig(x, eta_val=1, epsilon=0):
#    eta_val = eta(x,eta_0)
    Sideal = -R*(x * np.log(x) + (1 -x) * np.log(1-x))
    xpin = (x - epsilon) / (1 - epsilon)
    Sleft = R*(-xpin * np.log(xpin) - (1 - xpin) * np.log(1 - xpin)) *(1-epsilon)
    Sright = R*((eta_val) * np.log(eta_val) - (x) * np.log(x) - (eta_val - x) * np.log(eta_val - x))
    S = (Sleft * Sright) / Sideal
#    S = R*(-(xpin) * np.log(xpin) - (eta_val - xpin) * np.log(eta_val - xpin)) * (1- epsilon)
#    S = R*(-x * np.log(x) - (1 - x) * np.log(1 - x))
    return(S)

def dSconfig(S,x):
    dS = np.gradient(S) / np.gradient(x)
    return(dS)


fig, ((ax1, ax4),(ax2,ax5),(ax3,ax6)) = plt.subplots(3,2, sharex='col',gridspec_kw={'hspace': 0})
left_axes = [ax1,ax2,ax3]
right_axes = [ax4,ax5,ax6]

plt.subplots_adjust(left=0.50, bottom=0.10)
ax7 = fig.add_subplot(111)
axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7]
ax7.set_position(pos=[0.11,0.52,0.26,0.3])

# S1 = Sconfig(x, eta_1_init, epsilon_1_init) # Initialise with eta=1

voltage1, = ax1.plot(df_bw['x']*initial_cap, df_bw['VV_' + suffix], lw=2, label='model') # Initial plt
enthalpy1, = ax2.plot(df_bw['x']*initial_cap, -df_bw['dH_' + suffix]/F, lw=2, label= 'model')
entropy1, = ax3.plot(df_bw['x']*initial_cap, (T * df_bw['dS_' + suffix])/(F*1000), lw=2, label = 'model')

sublatt1, = ax4.plot(df_bw['x'], df_bw['n1'], lw=2,label ='Pores') # Initial plt
sublatt2, = ax4.plot(df_bw['x'], df_bw['n2'], lw=2,label ='Layers') # Initial plt
x1, = ax4.plot(df_bw['x'], df_bw['x'], lw=2,linestyle='--',color='grey',label='SS') # Initial plt
op1, = ax5.plot(df_bw['x'], df_bw['op'], lw=2,label ='Order param',linestyle='-') # Initial plt
S1, = ax6.plot(df_bw['x'], df_bw['S_'+suffix], lw=2,label ='Entropy',linestyle='-') # Initial plt
dQdV1, = ax7.plot(df_bw['VV_'+suffix], df_bw['dxdmu_'+suffix], lw=2,label ='dQ/dV',linestyle='-') # Initial plt

df_bw['op_null'] = 0
x = np.linspace(0,1,len(df_bw))
Sss =  Sconfig(x)
df_bw['S,ss'] = Sss
df_bw['dS,ss'] = dSconfig(Sss,x)

ax5.plot(df_bw['x'],df_bw['op_null'],linestyle='--',color='grey')
entropyss, = ax3.plot(df_bw['x']*initial_cap,(T * df_bw['dS,ss'])/(F*1000),linestyle='--',color='grey')
ax6.plot(df_bw['x'],df_bw['S,ss'],linestyle='--',color='grey')

ax3.set_xlabel('Capacity (mAh/g)')
ax6.set_xlabel('Sodiation fraction x')
ax7.set_xlabel('Voltage (V vs. Na)')

ax1.set_ylabel('OCV (V vs. Na)')
ax2.set_ylabel('-dH/dx (eV per Na)')
ax3.set_ylabel('T dS/dx (eV per Na)')

ax4.set_ylabel('Occupancy')
ax5.set_ylabel('Order parameter')
ax6.set_ylabel('S (J/mol/K)')

ax7.set_ylabel('dQ/dV (arb. units)')


for ax in left_axes:
    ax.autoscale()

for ax in right_axes:
    ax.yaxis.tick_right()
    ax.autoscale()
# ax.margins(x=0) # Whitespace padding

ax5.set_ylim([-1.3,1.3])
ax7.set_xlim([-0.2,0.3])

axcolor1 = 'lightgoldenrodyellow'
axcolor2 = 'darksalmon'

slider_spacing = 0.07
pad = 0.05
left = 0.10
width = 0.2
height = 0.015

axeta = plt.axes([left, 0.33+pad, width, height], facecolor=axcolor1) # Boxes for the slider. Args: left, bottom, width, height.
axE01 = plt.axes([left, 0.28+pad, width, height], facecolor=axcolor1) # Boxes for the 
axbeta = plt.axes([left, 0.23+pad, width, height], facecolor=axcolor2) # Boxes for the 
axalpha =plt.axes([left, 0.18+pad, width, height], facecolor=axcolor2) # Boxes for the
axdelE =plt.axes([left, 0.13+pad, width, height], facecolor=axcolor2) # Boxes for the
axL =plt.axes([left, 0.08+pad, width, height], facecolor=axcolor2) # Boxes for the
axg =plt.axes([left, 0.03+pad, width, height], facecolor=axcolor2) # Boxes for the 

# Calls slider itself. Args: ax, label, initial value, max value.
seta = Slider(axeta, 'Cap', 180.0, 500.0, valinit = initial_cap, valstep = 10.0)
setE01 = Slider(axE01, 'E0,1', -1.0, 0.2, valinit = -0.25, valstep = delta_eta)
setbeta = Slider(axbeta, 'beta', 0.25, 3.0, valinit = 2.0, valstep = delta_eta)
setalpha = Slider(axalpha, 'alpha', -2.5, -0.5, valinit = -1.0, valstep = delta_eta) 
setdelE =  Slider(axdelE, 'delE', -1.0, 0.0, valinit = -0.25, valstep = delta_eta)
setL =  Slider(axL, 'L', 0.10, 0.90, valinit = 0.25, valstep = delta_eta)
setg =  Slider(axg, 'g', -0.07, 0.02, valinit = -0.047, valstep = 0.002) 

def get_experimental_data():
    # accesses all the relevant plots.
    file_list = os.listdir(experimental_dir)
    entropy_file = [experimental_dir + f for f in file_list if f.endswith('entropy.csv')][0]    
    df_expt = pd.read_csv(entropy_file, encoding='latin')
    df_expt['Normalised cap [mAh/g]'] = df_expt['Charge/Discharge [mAh]']/mass
    return(df_expt)
    

def get_experimental_dqdv():
    dqdv = DQDV()
    dqdv.OCV_to_dQdV()
    ocv_df = dqdv.ocv_df_dis
    return(ocv_df)

def plot_experimental_data(df_expt):
    df_expt['M1 Enthalpy [J mol-1]'] = - F * df_expt['OCV [V]   '] * 1000 + T * df_expt['M1 Entropy [J mol-1 K-1]']    
    ax1.plot(df_expt['Charge/Discharge [mAh]']/mass,df_expt['OCV [V]   '] * 1000,linestyle='--',color='tomato',markersize = 2,label='Expt',marker='o') 
    ax2.plot(df_expt['Charge/Discharge [mAh]'][1:-2]/mass,-df_expt['M1 Enthalpy [J mol-1]'][1:-2]/96.485,linestyle='--',markersize=2,marker='o',label='Expt',color='tomato')
    ax3.plot(df_expt['Charge/Discharge [mAh]'][1:-2]/mass,(df_expt['M1 Entropy [J mol-1 K-1]'][1:-2] * T)/96.485,linestyle='--',markersize=2,marker='o',label='Expt',color='tomato')

df_expt = get_experimental_data()
plot_experimental_data(df_expt)
df_dqdv = get_experimental_dqdv()
ax7.plot(df_dqdv['OCV'],-df_dqdv['dQdV'],label='Expt')
max_cap = float(df_expt['Normalised cap [mAh/g]'].max()) # Maximum obtained capacity in mAh/g

for ax in [ax1,ax2,ax3,ax4,ax7]:
    ax.legend(fontsize=11)
    
#            ax1.plot(df_e['Charge/Discharge [mAh]'][1:-2]/mass,df_e['%s Entropy [J mol-1 K-1]'%method][1:-2]*288/96.485,linestyle=':',markersize=5,marker='^',color='grey',label='Entropy')    

# Auto updates values and y-axis scale.
def update(val):
    cap = seta.val
    beta = setbeta.val
    E01 = setE01.val
    alpha = setalpha.val
    delE = setdelE.val
    L = setL.val
    g = setg.val    
    arg_dict['E0'] = [E01]
    arg_dict['beta4'] = [beta]
    arg_dict['alpha4'] = [alpha]
    arg_dict['delE'] = [delE]
    arg_dict['L'] = [L]
    arg_dict['g2'] = [g]    
    two_layer = Two_layer(arg_dict) # Instantiate the class    
    df_bw = two_layer.solution() # Run all thermodynamics
    suffix = str(format(arg_dict[arg_dict['label']][0],'.3f'))    
        
#    S1 = Sconfig(x,eta)    
    voltage1.set_ydata(df_bw['VV_'+suffix])
    enthalpy1.set_ydata(-df_bw['dH_'+suffix]/F)
    entropy1.set_ydata((T * df_bw['dS_'+suffix])/(F*1000))
    voltage1.set_xdata(df_bw['x']*cap)
    enthalpy1.set_xdata(df_bw['x']*cap)
    entropy1.set_xdata(df_bw['x']*cap)
    entropyss.set_xdata(df_bw['x']*cap)    
    sublatt1.set_ydata(df_bw['n1']) # Initial plt
    sublatt2.set_ydata(df_bw['n2']) # Initial plt
    op1.set_ydata(df_bw['op']) # Initial plt
    S1.set_ydata(df_bw['S_'+suffix]) # Initial plt
    dQdV1.set_xdata(df_bw['VV_'+suffix])
    normalisation_factor = cap / max_cap
    dQdV1.set_ydata(df_bw['dxdmu_'+suffix])

#    freq = sfreq.val
#    l.set_ydata(amp*np.sin(2*np.pi*freq*t))
    fig.canvas.draw_idle()
#    for ax in axes:
#        ax.autoscale()


seta.on_changed(update)
setbeta.on_changed(update)
setE01.on_changed(update)
setalpha.on_changed(update)
setdelE.on_changed(update)
setL.on_changed(update)
setg.on_changed(update)


# resetax = plt.axes([0.8, 0.025, 0.1, 0.04]) # ?
# button = Button(resetax, 'Reset', color=axcolor1, hovercolor='0.975')

# def reset(event):
#    sfreq.reset()
#    samp.reset()
# button.on_clicked(reset)

#plt.tight_layout()
plt.show()

# if __name__ == '__main__':
