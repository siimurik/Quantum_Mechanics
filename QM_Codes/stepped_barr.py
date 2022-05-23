# Siim Erik Pugal
#===============================================================================
import numpy as np
import matplotlib.pyplot as plt

hbar = 1.055E-34
m   = 9.109E-31
E   = -0.8
U0  = -1.5
k1 = -np.sqrt(2*m*E)/hbar
k2 = -np.sqrt(2*m*(U0-E))/hbar

def P(x):
    P = 4*k1**2/(k1**2+k2**2)*np.exp(-2*k2*x)
    return P

x = np.linspace(-10, 10, 10000)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x, P(x))
#ax.set_title("Amdahl's law")
#ax.set_xlabel('Number of processors')
#ax.set_ylabel('Speedup')
plt.grid()
#plt.legend()
plt.show()
