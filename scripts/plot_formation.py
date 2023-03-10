import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

mpl.rcParams['lines.linewidth'] = 2.5

font = {'size': 18}

mpl.rc('font', **font)
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams.update({'errorbar.capsize': 2})
mpl.rcParams['legend.handlelength'] = 0

df = pd.read_csv('../Cathode/samsung_cathode_galvan_27mAg_00.txt',skiprows=12)
#df = pd.read_csv('../Cathode/samsung_cathode_00_13,5_refresh.txt',skiprows=12)
mass = df['I[A]'].max() / df['I[A/kg]'].max()

df['Cyc-Count-diff'] = df['Cyc-Count'].diff()
print('mass = ',mass)
df['Ah[Ah]'] = df['Ah[Ah/kg]'] * mass
df['Ah-Dis[Ah]'] = df['Ah-Dis[Ah/kg]'] * mass
df['Ah-Ch[Ah]'] = df['Ah-Ch[Ah/kg]'] * mass

last_index = 0
changes = df[df['Cyc-Count-diff'] != 0].index.tolist()
changes.append(df['Cyc-Count'].index[-2])


for cycle_no,index in enumerate(changes[1:]):
    last_index = changes[cycle_no]
    print(last_index,index)
    cycle_df = df.loc[last_index:index]

    if cycle_no == 1:
        max_initial_cap = cycle_df['Ah[Ah]'].max()
        min_initial_cap = cycle_df['Ah[Ah]'].min()
#        max_initial_cap = cycle_df['Ah[Ah]'].max()
#        min_initial_cap = cycle_df['Ah[Ah]'].min()                
        correction = min_initial_cap
        print('max_initial_cap=',max_initial_cap)
    elif cycle_no > 1:
#        min_cycle_cap = cycle_df['Ah[Ah/kg]'].min()
        min_cycle_cap = cycle_df['Ah[Ah]'].min()        
        correction = min_cycle_cap
#        correction = 0
#        print('max_cycle_cap=',max_cycle_cap)
    else:
        correction = 0
    
#    cycle_df = cycle_df[cycle_df['I[A/kg]'] != 0] # Exclude OCV
    cycle_df['Capacity'] = cycle_df['Ah[Ah]'] - correction
    cycle_df_c = cycle_df[(cycle_df['I[A/kg]']) > 9]
   
#    plt.plot(-cycle_df_c['Capacity'],cycle_df_c['U[V]'],label='Charge,'+str(cycle_no+1))

    
    cycle_df_d = cycle_df[((cycle_df['I[A/kg]']) < -9) | ((cycle_df['I[A/kg]']) > 9)]
#    plt.plot(cycle_df_d['~Time[h]'],cycle_df_d['I[A/kg]'],label='Charge,'+str(cycle_no+1),linestyle='',marker='o')
    if cycle_no in [1,2,3]:
        line_list = ['-','dashed','dotted','','--']
        colours=['royalblue','grey','black','','salmon']
        plt.plot(cycle_df_d['Capacity']*1000,cycle_df_d['U[V]'],label='Cycle '+str(cycle_no),linestyle=line_list[cycle_no-1],color=colours[cycle_no-1])
        cycle_df_d.to_csv('cathode_galvanostatic_cycle_{}.csv'.format(cycle_no))
    last_index = index

plt.xlabel('Capacity (mAh)')
plt.ylabel('Voltage (V) vs. Li') 
plt.legend(handlelength=2)
plt.savefig('3_cycles_cathode.png')
plt.show()
'''
print(changes)
plot, = plt.plot(-df['Ah[Ah/kg]'].iloc[52491:],df['U[V]'].iloc[52491:])
plt.legend()
plt.show()
'''
