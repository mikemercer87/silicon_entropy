import pandas as pd
import matplotlib.pyplot as plt

x = [100,80,60,15,10]
y = [0.446,0.418,0.399,0.215,0.213]

plt.plot(x,y,marker='x',linestyle='--')
plt.xlabel('SiOx mass fraction')
plt.ylabel('Approximate peak voltage (V)')
plt.show()