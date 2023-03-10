import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter

'''
Script that handles postprocessing to get dQ/dV files.
'''

class DQDV():
    def __init__(self,input_file = None):
        self.dqdv = pd.DataFrame()
        self.entropy_path = '../entropy_data/'
        self.figures_path = '../figures/'
        self.galvan_path = '../'
        self.input_file = input_file
      
    def data_reduction(self, input_df, slice_window):
        # Cuts down the number of data points prior to smoothing.
        # Generates new dataframe, uses high resolution OCV channels.
        dqdv = pd.DataFrame()
        dqdv['SOC'] = self.log_rolling_average(input_df['SOC'],slice_window)
        self.ocv_key = str([k for k in input_df if k.startswith('OCV')][0])
        print('ocv_key='+self.ocv_key)
        dqdv['U[V]'] = self.log_rolling_average(input_df[self.ocv_key],slice_window)
        dqdv['Time'] = input_df['~Time[h]']
        dqdv['I[A]'] = input_df['I[A]']
        return(dqdv)
        
    def log_rolling_average(self, df_col, slice_window):
        rolling_avg = df_col.rolling(int(slice_window),center=True).mean().dropna()
        print((df_col.rolling(int(slice_window),center=True)), 'ROLLING')
        return(rolling_avg)
            
    def differentiate(self, df, ir = False, state='dis'):
        # Performs differentiation of the dataset.
        dx = float(df['SOC'].iloc[1]) - float(df['SOC'].iloc[0])
        dV = float(df['U[V]'].iloc[1]) - float(df['U[V]'].iloc[0])

        deriv = np.gradient(df['U[V]'],edge_order=2)
        deriv2 = np.gradient(df['SOC'],edge_order=2)             

        df['dE/dx'] = -deriv / dx
        df['dx/dE'] = 1 / (df['dE/dx'])
        return(df['dx/dE'])

    def read_files(self, key='Galvan_C50', cells = 'all'):
        # Export processed preformation into csv. If cells is entered as a list, only those ones are processed.
        # Plots U as a function of time, for all cycles.
        self.df = pd.read_csv(self.input_file, skiprows=12, sep=',')         
        print(list(self.df))
        self.df_split()

    def df_split(self):
        self.all_lines = list(self.df['Line'].drop_duplicates().values)
        self.all_lines.reverse()
        print('all_lines=', self.all_lines)
        self.split_dfs = {}
        self.full_discharge = self.df[self.df['Line'] == 3] # You would need to customise based on your test plan.
        self.full_charge = self.df[self.df['Line'] == 5]
        self.cccv = self.df[self.df['Line'] == 4]
        self.max_dc = self.full_discharge['Ah-Dis[Ah/kg]'].max()
        self.max_ch = self.full_charge['Ah-Ch[Ah/kg]'].max()
        self.min_dc = self.full_discharge['Ah-Dis[Ah/kg]'].min()
        self.min_ch = self.full_charge['Ah-Ch[Ah/kg]'].min()       
        self.cccv_corr = self.cccv['Ah-Dis[Ah/kg]'].max() - self.max_dc
        print('max_ch =', self.max_ch, 'max_dc=', self.max_dc)
        print('correction=', self.cccv_corr)
        
        for line in self.all_lines:
            self.split_dfs[line] = self.df[self.df['Line'] == line]
            
    def plot(self,kind='lithiation',color='blue',label='text'):
        print('kind=', kind)
        linestyles = ['-','dashed','dashed',':',':']
        colors=['royalblue','grey','salmon','olive','purple']

        for window in [300]: # Window is the number of data points to use in data reduction.
            last_ch = 0
            last_dc = 0
            self.read_files()            
#            for index, count in enumerate(range(1,self.df['Count'].max()+1)):
            for index, count in enumerate(range(1,2)):      # Edit to choose a different cycle.          
                self.df_cyc = self.df[self.df['Count'] == count]

                self.df_cyc['SOC'] = ((-self.df_cyc['Ah[Ah/kg]'] + self.df_cyc['Ah[Ah/kg]'].iloc[0]) / (self.max_dc - self.df_cyc['Ah[Ah/kg]'].iloc[0]))
#                self.df_cyc['SOC'] = (self.df_cyc['Ah[Ah/kg]']) / (self.max_dc)
                
                if kind == 'lithiation':
                    self.df_cyc = self.df_cyc[self.df_cyc['U[V]'] < 1.0] # Exclude constant voltage
                    self.df_cyc = self.df_cyc[self.df_cyc['I[A]'] < 1e-9] # Discharge cycles only.                    
                else:
                    self.df_cyc = self.df_cyc[self.df_cyc['U[V]'] > 0.020] # Exclude constant voltage
#                    self.df_cyc = self.df_cyc[self.df_cyc['I[A]'] ] # Charge cycles only.                                        

                self.dqdv = self.data_reduction(self.df_cyc, slice_window = window)
                self.dqdv['dqdv'] = self.differentiate(self.dqdv)
                if kind == 'lithiation':                    
                    self.dqdv_red =self.dqdv[self.dqdv['dqdv'] > 0]
                else:                    
                    self.dqdv_red =self.dqdv[self.dqdv['dqdv'] < 0]                            
            
                plt.ylim([-100,100])
                plt.xlabel('SOC')

                if kind == 'lithiation':
#                    plt.plot(self.dqdv_red['U[V]'].iloc[1:-1], -self.dqdv_red['dqdv'].iloc[1:-1], label='cycle='+str(count)+', ' + kind,color=colors[index+1],linestyle=linestyles[index+1])
                    plt.plot(self.dqdv_red['U[V]'].iloc[1:-1], -self.dqdv_red['dqdv'].iloc[1:-1], label=label+',l',color=colors[index],linestyle=linestyles[index+1])
                else:
                    plt.plot(self.dqdv_red['U[V]'].iloc[1:-1], -self.dqdv_red['dqdv'].iloc[1:-1], label=label+',d',color=color,linestyle=linestyles[index])
                             
                plt.ylabel('dQ/dV (arb units)')
                plt.xlabel('Voltage (V) vs. Li')
                
#                self.fig2.savefig('dVdQ_%s_ch.png' % self.cell_number)
    def OCV_to_dQdV(self):
        df_d = pd.read_csv(self.entropy_path + 'entropydischarge_hc_00entropy.csv')
        df_d['OCV'] = df_d['OCV [V]   '] * 1000 # Converts to volts
        df_d['x'] = df_d['Charge/Discharge [mAh]'] / df_d['Charge/Discharge [mAh]'].iloc[-1]
        df_d['dQdV'] = np.gradient(df_d['x'])/np.gradient(df_d['OCV'])
#        plt.plot(df_d['OCV'],df_d['dQdV'],linestyle='',marker='o',label='OCV, sodiation')
        self.ocv_df_dis = df_d # OCV from discharge
        df_c = pd.read_csv(self.entropy_path + 'entropycharge_hc_00entropy.csv')
        df_c['OCV'] = df_c['OCV [V]   '] * 1000 # Converts to volts
        df_c['x'] = df_c['Charge/Discharge [mAh]'] / df_c['Charge/Discharge [mAh]'].iloc[-1]
        df_c['dQdV'] = -np.gradient(df_c['x'])/np.gradient(df_c['OCV'])
        self.ocv_df_ch = df_c
        
#        plt.plot(df['OCV'],df['dQdV'],linestyle='',marker='^',label='OCV, desodiation')        
                
if __name__ == '__main__':
    dqdv = DQDV(input_file = '../amte_grsiox_12,3mAg_03.txt')
    dqdv.plot(kind='lithiation',color='orange',label='SiOx')
    dqdv.plot(kind='delithiation',color='dodgerblue',label='SiOx')

    dqdv = DQDV(input_file = '../amte_gr_12,3mAg_04.txt')
    dqdv.plot(kind='lithiation',color='grey',label='Gr')
    dqdv.plot(kind='delithiation',color='darkblue',label='Gr')    
    
#    dqdv.OCV_to_dQdV()
    plt.legend()
    plt.savefig('plot_overlaid.png',dpi=500)
    plt.show()    

    plt.clf()


    
