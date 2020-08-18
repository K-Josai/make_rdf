import pandas as pd
import matplotlib.pyplot as plt
# import numpy as np


df = pd.read_csv("g_of_r.txt", delim_whitespace=True, header=None)
plt.ylim(0,10)
plt.xlim(0,5)
plt.plot(df[0], df[1])
plt.show()
