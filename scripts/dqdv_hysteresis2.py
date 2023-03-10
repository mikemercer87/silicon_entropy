import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.interpolate import interp1d

'''
Script that handles postprocessing to get dQ/dV files.
'''

class DQDV():
    def __init__(self):
#        Plotter.__init__(self) # Inherit variables from Plotter.
#        self.plotter_inst = Plotter() # Instantiate this class.
        self.dqdv = pd.DataFrame()
      
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
        dqdv['U[V]'] = self.log_rolling_average(input_df['OCV'],slice_window)
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
        self.df = pd.read_csv('../SiOx/amte_grsiox_hysteresis_00.txt',skiprows=12,sep=',')
        
        print(list(self.df))

        self.df_split()

    def df_split(self):
        self.all_lines = list(self.df['Line'].drop_duplicates().values)
#        self.all_lines.reverse()
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
            

    def plot(self):
        self.read_files()
#        self.all_lines.reverse()        
        for window in [25]:
            last_ch = 0
            last_dc = 0
            dash = '--'
            print(self.all_lines)
    
            for line in self.all_lines[1:]:
#            for line in [7,10,13,16]:
#                print('line=',line)
 #               if line == 4 or line == 5:
#                colour = 'r'
#                dash = 'dotted'
                label='None'
                if line == 6 or line == 7:
                    colour = 'red'
                    dash = '--'
                    label = '1'
                elif line == 9 or line == 10:
                    colour ='skyblue'
                    dash = '--'
                    label = '2'
                elif line == 12 or line == 13:
                    colour = 'grey'
                    dash = '--'
                    label='3'
                elif line == 15 or line == 16:
                    colour = 'black'
                    dash = '--'
                    label='4'
#                elif line == 17 or line == 18:
#                    colour = 'mediumseagreen'
#                    dash = '-'                    
#                elif line == 20 or line == 21:
#                    colour = 'green'
#                    dash = '-'
#                elif line == 23 or line == 24:
#                    colour = 'darkolivegreen'
#                    dash = '-'                        
                else:
                    colour = 'royalblue'
                    dash = 'None'                    

#                colour = 'r'
#                dash = '--'

                
                self.df = self.split_dfs[line]
                self.df['OCV']=self.df['OCV01[V]']
                if self.df['I[A]'].iloc[2] > 0:  # Charge
                    self.df['SOC'] = (-self.df['Ah[Ah/kg]']) / (self.max_ch) 
#                    if line == 5 or line == 13: # Full charge cycle:
#                        self.df['SOC'] = self.df['SOC'] - self.cccv_corr / self.max_ch # Correction for extra capacity gained during cccv   
#                    last_ch = self.df['SOC'].max()
                else:
                    self.df['SOC'] = -self.df['Ah[Ah/kg]'] / self.max_ch                    
#                    self.df['SOC'] = (-self.df['Ah[Ah/kg]'] + self.df['Ah[Ah/kg]'].iloc[0]) / (self.max_dc - self.df['Ah[Ah/kg]'].iloc[0])

#                if self.df['SOC'].max() > 1:
#                    self.df['SOC'] = self.df['SOC'] - (self.df['SOC'].max() - 1)

#                if self.df['SOC'].min() > 0:
#                    self.df['SOC'] = self.df['SOC'] - (self.df['SOC'].min())

#                if line > 3:
#                    self.df['SOC'] = self.df['SOC'] - self.cccv_corr
#                    last_dc = self.df['SOC'].min()
#                    kwindow = window /4000
                self.dqdv = self.data_reduction(self.df, slice_window = window)
                self.dqdv['dqdv'] = self.differentiate(self.dqdv)
                self.dqdv_red =self.dqdv[self.dqdv['dqdv'] > 0]
                if line == 3:
                    self.dqdv_red.to_csv('lithiation_25.csv')
                elif line == 5:
                    self.dqdv_red.to_csv('delithiation_25.csv')
#                plt.xlim([0.02,0.25])
#                plt.xlim([0.10,0.24])                
#                plt.ylim([-10000,10000])
                plt.xlabel('SOC')
                if label != 'None':
                    if self.df['I[A]'].iloc[5] > 0:  # Charge
                        plt.plot(self.dqdv_red['U[V]'], self.dqdv_red['SOC'], label='line='+label,color=colour,linestyle=dash)
                    else:
                        plt.plot(self.dqdv_red['U[V]'], self.dqdv_red['SOC'], label='line='+label,color=colour,linestyle=dash)
                plt.ylabel('SOC')
        plt.legend()
        plt.savefig('plot_04_25.png',dpi=500)
        plt.show()    

        plt.clf()

#                self.fig2.savefig('dVdQ_%s_ch.png' % self.cell_number)                
                
if __name__ == '__main__':
    dqdv = DQDV()
    dqdv.plot()
