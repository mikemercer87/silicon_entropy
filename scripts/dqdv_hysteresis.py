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
    def __init__(self,input_file):
#        Plotter.__init__(self) # Inherit variables from Plotter.
#        self.plotter_inst = Plotter() # Instantiate this class.
        self.dqdv = pd.DataFrame()
        self.entropy_path = '../entropy_data/'
        self.figures_path = '../figures/'
        self.galvan_path = '../galvanostatic/'
        self.input_file = input_file
      
    def cubic_smoother(self, df, option = 'SOC', n_points = 5000):
        # Returns two smoothed and interpolated grids, x and the variable. Npoints is the number of points to smooth over.
        df_smoothed = pd.DataFrame() # Empty dataframe for the smoothed data.
        
        if option == 'SOC':
            smooth= smooth_fn(new_x)
            df_smoothed['SOC_%d' % n_points] = new_x
            df_smoothed['dx/dE_%d' % n_points] = smooth

        return(df_smoothed)

    def data_reduction(self, input_df, slice_window):
        # Cuts down the number of data points prior to smoothing.
        # Generates new dataframe, uses high resolution OCV channels.
        dqdv = pd.DataFrame()
        dqdv['SOC'] = self.log_rolling_average(input_df['SOC'],slice_window)
        self.ocv_key = str([k for k in input_df if k.startswith('OCV')][0])
        print('ocv_key='+self.ocv_key)
        dqdv['U[V]'] = self.log_rolling_average(input_df[self.ocv_key],slice_window)
        dqdv['Time']=input_df['~Time[h]']
        dqdv['I[A]']=input_df['I[A]']
        return(dqdv)

    def interpolation(self,interp_percent = 2):
        self.dqdv_dis_s = pd.DataFrame() # smoothed results
        self.dqdv_ch_s = pd.DataFrame()
        self.smooth_dfs = [self.dqdv_dis_s,self.dqdv_ch_s]  
        
        for df, smooth_df in zip(self.raw_dfs,self.smooth_dfs):
            df_sorted = df.sort_values(by = 'SOC')
            min_SOC = df_sorted['SOC'].min()
            max_SOC = df_sorted['SOC'].max()
            min_E = df_sorted['U[V]'].min()
            max_E = df_sorted['U[V]'].max()            
            n_points = len(df_sorted['SOC'])
            new_x = np.linspace(start = min_SOC, stop = max_SOC, num = (interp_percent/100) * n_points)
            new_E = np.linspace(start = min_E, stop = max_E, num = (interp_percent/100) * n_points)            
            soc_array = df_sorted['SOC'].values
            dqdv_array = df_sorted['dx/dE'].values
            E_array = df_sorted['U[V]'].values
            smooth_dqdv = interp1d(x = soc_array.tolist(), y = dqdv_array.tolist(), kind='linear') # Try only single interpolation!
#            smooth_dqdv = interp1d(x = soc_array.tolist(), y = dqdv_array.tolist(), kind='linear')            
            smooth_df['SOC'] = new_x
            smooth_df['dx/dE'] = smooth_dqdv(new_x)
            smooth_df['U[V]'] = E_array
        
    def log_rolling_average(self, df_col, slice_window):
        # Amended to use savitsky-golay
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
        self.df = pd.read_csv(self.input_file,skiprows=12,sep=',')
#        self.df = pd.read_csv('galvanostatic/HC_formation_00.txt',skiprows=12,sep=',')        
        
        print(list(self.df))

        self.df_split()

    def df_split(self):
        self.all_lines = list(self.df['Line'].drop_duplicates().values)
        self.all_lines.reverse()
        print('all_lines=', self.all_lines)
        self.split_dfs = {}
        self.full_discharge = self.df[self.df['Line'] == 3]
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
        linestyles = ['-','dashed','dashed',':',':']
        colors=['royalblue','grey','salmon','olive','purple']
        self.read_files()
        for window in [300]:
            last_ch = 0
            last_dc = 0
            dash = '--'
            
            for index, count in enumerate(range(1,self.df['Count'].max()+1)):
#            for index, count in enumerate(range(1,2)):                
                self.df_cyc = self.df[self.df['Count'] == count]

                self.df_cyc['SOC'] = ((-self.df_cyc['Ah[Ah/kg]'] + self.df_cyc['Ah[Ah/kg]'].iloc[0]) / (self.max_dc - self.df_cyc['Ah[Ah/kg]'].iloc[0]))*2
#                self.df_cyc['SOC'] = (self.df_cyc['Ah[Ah/kg]']) / (self.max_dc)
                
                if kind == 'lithiation':
                    self.df_cyc = self.df_cyc[self.df_cyc['U[V]'] < 1.0] # Exclude constant voltage
                    self.df_cyc = self.df_cyc[self.df_cyc['I[A]'] < -1e-9] # Discharge cycles only.                    
                else:
                    self.df_cyc = self.df_cyc[self.df_cyc['U[V]'] > 0.032] # Exclude constant voltage
#                    self.df_cyc = self.df_cyc[self.df_cyc['I[A]'] > -1e-3] # Charge cycles only.                                        

                self.dqdv = self.data_reduction(self.df_cyc, slice_window = window)
                self.dqdv['dqdv'] = self.differentiate(self.dqdv)
                if kind == 'lithiation':                    
                    self.dqdv_red =self.dqdv[self.dqdv['dqdv'] > 0]
                else:                    
                    self.dqdv_red =self.dqdv[self.dqdv['dqdv'] < 0]                            

#                plt.xlim([0.02,0.25])
#                plt.xlim([0.10,0.24])                
#                plt.ylim([-100,100])
                plt.xlabel('SOC')

#                plt.plot(self.df_cyc['~Time[h]'].iloc[:-1], self.df_cyc['OCV01[V]'].iloc[:-1], label='cycle='+str(count)+', ' + kind)
#                plt.plot(self.df_cyc['~Time[h]'].iloc[:-1], self.df_cyc['I[A]'].iloc[:-1]*1e6, label='cycle='+str(count)+', ' + kind)
                if kind == 'lithiation':
#                    plt.plot(self.dqdv_red['U[V]'].iloc[1:-1], -self.dqdv_red['dqdv'].iloc[1:-1], label='cycle='+str(count)+', ' + kind,color=colors[index],linestyle=linestyles[index])
                    plt.plot(self.dqdv_red['U[V]'].iloc[1:-1], -self.dqdv_red['SOC'].iloc[1:-1], label=str(count)+',l',color=color,linestyle=linestyles[index])
                else:
                    plt.plot(self.dqdv_red['U[V]'].iloc[1:-1], -self.dqdv_red['SOC'].iloc[1:-1], label=str(count)+',d',color=color,linestyle=linestyles[index+1])
                             
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
#    dqdv = DQDV(input_file = '../SiOx/amte_grsiox_27mAggalvan_01.txt')
#    dqdv.plot(kind='lithiation')
#    dqdv.plot(kind='delithiation')
    
#    dqdv = DQDV(input_file = '../SiOx/amte_grsiox_18,2mAg_01.txt')

#    dqdv = DQDV(input_file = '../SiOx/amte_grsiox_12,3mAg_03_partial.txt')
#    dqdv = DQDV(input_file = '../Gr/amte_gr_27mAg_05.txt')        
#    dqdv.plot(kind='lithiation',color='r',label='27 mA/g')
#    dqdv.plot(kind='delithiation',color='b',label='27 mA/g')


#    dqdv = DQDV(input_file = '../SiOx/amte_grsiox_12,3mAg_03.txt')
#    dqdv.plot(kind='lithiation',color='orange',label='SiOx')
#    dqdv.plot(kind='delithiation',color='dodgerblue',label='SiOx')

    dqdv = DQDV(input_file = '../SiOx/amte_grsiox_hysteresis_partial.txt')
    dqdv.plot(kind='lithiation',color='grey',label='Gr')
#    dqdv.plot(kind='delithiation',color='darkblue',label='Gr')    

    
 #   dqdv.OCV_to_dQdV()
    plt.legend()
    plt.savefig('plot_overlaid.png',dpi=500)
    plt.show()    

    plt.clf()


    
