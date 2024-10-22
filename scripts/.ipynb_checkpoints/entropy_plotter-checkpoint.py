import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import string
# from dqdv_proc import DQDV
from scipy.integrate import cumtrapz
# Directories.

def plot_formatting():
    # Reconfigure the plot display as you like (refer to Matplotlib documentation for details).
    mpl.rcParams['lines.linewidth'] = 1.5
    mpl.rcParams['lines.markersize'] = 4
    font = {'size': 18}
    mpl.rc('font', **font)
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams.update({'errorbar.capsize': 2})
    mpl.rcParams['legend.handlelength'] = 2
#********* change the path below to plot however many entropy plots you like.
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

entropy_path = '../data/' # You should save the csvs produced by Matlab in this folder.
figures_path = '../figures/'   # relative path where figures are saved
galvan_path = '../galvanostatic/' # this should contain your galvanostatic data. Keep the data file originally used to process and the .csv from the Matlab scripts with their default names.
temp_dirs = [entropy_path + 'AMTE_Gr_10SiOx/', entropy_path + 'AMTE_Gr_15SiOx/', entropy_path + 'AMTE_Gr_60SiOx/',entropy_path + 'AMTE_Gr_80SiOx/',entropy_path + 'AMTE_SiOx_SEI_Entropy/']

method = ['M1','M2','M3','M4','Bestfit']
#raw_keys = ['T','Smin','Smax','Hmin','Hmax','E1','E2','Estep']
#amplitudes_raw = {k: [] for k in raw_keys}
#processed_keys = ['T','Smin_avg','Smin_dev','Smax_avg','Smax_dev','Hmin_avg','Hmin_dev','Hmax_avg','Hmax_dev','E1_avg','E1_dev','E2_avg','E2_dev','Estep_avg','Estep_dev']
#amplitudes_proc = {k : [] for k in processed_keys} 

fig, (ax1, ax2, ax3) = plt.subplots(3,1, sharex = 'col', figsize=(6,14)) 

ax1.set_ylabel(r'${\Delta}H$' +' / kJ mol'+'$^{-1}$' +'K'+'$^{-1}$')
ax2.set_ylabel(r'${\Delta}S$' +' / J mol'+'$^{-1}$' +'K'+'$^{-1}$')
ax3.set_xlabel('OCV / V vs. Li')
ax3.set_ylabel('OCV / V vs. Li')
ax1.tick_params(which='both', length=10, width=1.5, direction='in',top=True,right=True)
ax1.tick_params(which='minor', length=5, width=1.0, direction='in',top=True,right=True,left=True,bottom=True)
ax2.tick_params(which='both', length=10, width=1.5, direction='in',top=True,right=True)
ax2.tick_params(which='minor', length=5, width=1.0, direction='in',top=True,right=True,left=True,bottom=True)
ax3.tick_params(which='both', length=10, width=1.5, direction='in',top=True,right=True)
ax3.tick_params(which='minor', length=5, width=1.0, direction='in',top=True,right=True,left=True,bottom=True)
# Physical constants.

left, bottom, width, height = [0.32, 0.18, 0.46, 0.3] #

F = 96485.3329
T_K = 298 # T in Kelvin. Note: you should customise if you change the temperature.
T = str(T_K - 273) # Temperature in degrees C, for labelling plots.

lower_df = pd.DataFrame()
upper_df = pd.DataFrame()

# plotting.


#def plot_wrt_cap():

#def plot_wrt_ocv():

for index,(directory) in enumerate(temp_dirs):
    file_list = os.listdir(directory)
    entropy_file_dict= {f.split('_')[1]: directory + f for f in file_list if f.endswith('entropy.csv')}
    basytec_file_dict = {f.split('_')[1]: directory + f for f in file_list if f.endswith('txt')}  # Puts all the files into a dictionary.
    print(list(entropy_file_dict))
    labels = ['10 %','15 %','60 %','80 %', 'SiOx'] # You can customise the legend labels here.
    for k in (entropy_file_dict):
        print(entropy_file_dict[k])
        entropy_file = entropy_file_dict[k]
        basytec_file = basytec_file_dict[k]
        print(entropy_file)
        df_e = pd.read_csv(entropy_file, encoding='latin')
        df_b = pd.read_csv(basytec_file, skiprows=12, encoding='latin')
        mass = df_b['Ah[Ah]'].iloc[-1] / df_b['Ah[Ah/kg]'].iloc[-1]  # Mass in kg.
        
        df_e['SOC'] = df_e['Charge/Discharge [mAh]'] / (df_e['Charge/Discharge [mAh]'].iloc[-1]) 

        df_e['M1 Enthalpy [J mol-1]'] = - F * df_e['OCV [V]   '] + T_K * df_e['M1 Entropy [J mol-1 K-1]']
        df_e['M2 Enthalpy [J mol-1]'] = - F * df_e['OCV [V]   '] + T_K * df_e['M2 Entropy [J mol-1 K-1]']
        df_e['M3 Enthalpy [J mol-1]'] = - F * df_e['OCV [V]   '] + T_K * df_e['M3 Entropy [J mol-1 K-1]']
        df_e['M4 Enthalpy [J mol-1]'] = - F * df_e['OCV [V]   '] + T_K * df_e['M4 Entropy [J mol-1 K-1]']

        df_e['M3 Enthalpy_Lower'] = df_e['M3 Entropy_Lower [J mol-1 K-1]'] * T_K
        df_e['M3 Enthalpy_Upper'] = df_e['M3 Entropy_Upper [J mol-1 K-1]'] * T_K
        

        for method in ['M1']: # Customise if you want a different fiiting methos.
            print('index =', index)
            cap_max = df_e['Charge/Discharge [mAh]'].max()
#            if index == 1:
#                df_e['Charge/Discharge [mAh]'] = df_e['Charge/Discharge [mAh]'] + cap_max
#            ax1.plot((df_e['Charge/Discharge [mAh]'][1:-1]) / (mass * 1000),df_e['%s Enthalpy [J mol-1]' % method][1:-1],linestyle='--',marker='o',color = colors[index], label = labels[index])
#            ax2.plot((df_e['Charge/Discharge [mAh]'][1:-1]) / (mass * 1000),df_e['%s Entropy [J mol-1 K-1]' % method][1:-1] * 1000,linestyle='--',marker='^',color = colors[index], label = labels[index])
#            ax3.plot((df_e['Charge/Discharge [mAh]'].iloc[1:-1]) / (mass * 1000),df_e['OCV [V]   '].iloc[1:-1] * 1000,linestyle='--',marker='o',color = colors[index],label = labels[index])
            ax1.plot(df_e['OCV [V]   '].iloc[1:-1]*1000,df_e['%s Enthalpy [J mol-1]' % method][1:-1],linestyle='--',marker='o',color=colors[index], label = labels[index])
            ax3.plot(df_e['OCV [V]   '].iloc[1:-1]*1000,df_e['OCV [V]   '].iloc[1:-1]*1000,linestyle='--',marker='o',color=colors[index],label=labels[index])
            ax2.plot(df_e['OCV [V]   '].iloc[1:-1]*1000,df_e['%s Entropy [J mol-1 K-1]' % method][1:-1]*1000,linestyle='--',marker='^',color=colors[index], label = labels[index])

        new_df = pd.DataFrame()
        new_df['OCV'] = df_e['OCV [V]   '] * 1000        
        new_df['Cap (mAh/g)'] = df_e['Charge/Discharge [mAh]']/mass
        new_df['deltaH (kJ/mol)'] = df_e['%s Enthalpy [J mol-1]'%method] *1000
        new_df['deltaS (J/mol/K)'] = df_e['%s Entropy [J mol-1 K-1]'%method]*1000
        new_df.to_csv('output.csv') # You can rename this manually.

        df_e.to_csv('data_entropy_%s_cap.csv' % index)           

ax3.legend()

if __name__ == '__main__':
    plot_formatting()
    plt.tight_layout()
    plt.savefig(figures_path + 'Enthalpy_entropy_ocv_anode_charge_cap.png') # Save the plot. Note: this will currently overwrite the existing plot. You would need to change the output file name randomly, or set up an automated naming system
    plt.show()

