import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

x       = np.arange(0,1,0.0001)
h       = 6.626068E-34 # m^2 kg / s
hbar    = h/(2*np.pi)
m       = 9.1093837E-31 # mass of electron
eV      = 1.60217663E-19 # 1 eV. Just in case
L_input = float(input("\nEnter barrier width in nm  = "))
L       = L_input*10**(-9)
U0_input= float(input( "Enter barrier height in eV = "))
U0      = U0_input*eV
#a=float(input("Enter efective barrier width = "))
ad      = L*np.sqrt(2*m*U0)/hbar
beta    = ad

def func(x):
    return np.tan(beta*np.sqrt(x))-np.sqrt(x*(1-x))/(x-0.5)

nmax    = int(np.trunc(beta/np.pi+0.5))+1
xdata   = np.zeros(nmax)

print("\n+----------------+----------------------+--------------------+",
      "\n|  Energy level  |  Infinite well (eV)  |  Finite well (eV)  |",
      "\n+----------------+----------------------+--------------------+")

xtest   = 0
for i in range(nmax):
    xdata[i]=(np.pi*(i+0.5)/beta)**2
    if (i==0):
        continue
    if (xdata[i-1]<0.5 and xdata[i]>0.5):
        xdata[i] = 0.5
        xtest    = 1
        continue
    if (xtest == 1):
        xdata[i]=(np.pi*(i-0.5)/beta)**2

lp  = 1.0001
rp  = 0.9999

for i in range(nmax):
    E_inf = (np.pi*hbar)**2/(2*m*L**2)*(1+i)**2/eV
    tout  = 0
    lv    = xdata[i]*lp
    if (i == nmax-1):
        rv=1.
    else:
        rv=xdata[i+1]*rp
    if (func(lv)*func(rv)>0):
        tout=1
    else:
        root=optimize.brentq(func,lv,rv,xtol=1.e-7,rtol=1.e-8)
    if (tout != 1):
        #print ("For energy level",i,"energy is:",root)
        print("\t", i+1, "\t   ", E_inf, "\t  ", root*U0/eV)
        #print("\t", i, "\t   ", np.round(E_inf, 6), "\t  ", np.round(root*U0/eV, 6))
