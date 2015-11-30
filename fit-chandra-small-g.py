tab=[[0.6, "+0.061", "+0.1827", 3.65], [0.8, -0.113, "+0.2079", 4.16], [1.0, -0.214, "+0.2138", 4.28], [1.2, -0.252, "+0.2047", 4.09], [1.4, -0.243, "+0.1862", 3.72], [1.6, -0.207, "+0.1634", 3.27], [1.8, -0.16, "+0.1404", 2.81], [2.0, -0.114, "+0.1192", 2.38], [3.0, "+0.027", "+0.0555", 1.11], [4.0, "+0.073", "+0.0313", 0.63], [5.0, "+0.092", "+0.0200", 0.4], [20.0, "+0.125", "+0.0013", 0.03], [30.0, "+0.125", 0.0006, 0.01], [40.0, "+0.125", 0.0003, "6e-3"], [50.0, "+0.125", 0.0002, "4e-3"], [100.0, "+0.125", "5e-5", "1e-3"], ["\\infty", "+0.125", "#ERROR", "#ERROR"]]
import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

a = np.array(tab[:-1], dtype=float)
x = a[:, 0]
g = a[:, 1]
z = np.polyfit(np.exp(-x), g, 2)
z[-1] = 0.125
p = np.poly1d(z)
xx = np.linspace(0.0, 10.0, 200)
fig, ax = plt.subplots(1, 1)
ax.plot(x, g, 'o')
ax.plot(xx, p(np.exp(-xx)), '-')
ax.set_xlim(0.0, 10)
ax.set_ylim(-0.5, 0.5)
ax.set_title(','.join(['{:.3f}'.format(_) for _ in z]))
figname = sys.argv[0].replace('.py', '.pdf')
fig.set_size_inches(4, 4)
fig.tight_layout()
fig.savefig(figname)
print(figname)
