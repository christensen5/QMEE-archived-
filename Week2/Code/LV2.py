""" The typical Lotka-Volterra Model simulated using scipy. This version differs from LV1 in that it accepts
    LV parameters as user-specified arguments, and saves the resulting plot to disk."""
__author__ = 'Alexander Christensen'
__version__ = '0.0.1'

import sys
import scipy as sc 
import scipy.integrate as integrate
import pylab as p  # Contains matplotlib for plotting
from matplotlib.offsetbox import AnchoredText

# import matplotlip.pylab as p #Some people might need to do this

def dCR_dt(pops, t=0):
    """ Returns the growth rate of predator and prey populations at any 
    given time step """
    
    R = pops[0]
    C = pops[1]
    dRdt = r*R*(1-(R/K)) - a*R*C
    dydt = -z*C + e*a*R*C
    
    return sc.array([dRdt, dydt])

""" The remaining part of the function carries out the integration step in simulating LV populations, and plots the
    result. Takes user-supplied inputs for parameters:
    r = resource growth rate
    a = consumer resource rate
    z = consumer mortality rate
    e consumer production efficiency"""

# Store parameters from user input:
r = float(sys.argv[1]) # Resource growth rate
a = float(sys.argv[2]) # Consumer search rate (determines consumption rate)
z = float(sys.argv[3]) # Consumer mortality rate
e = float(sys.argv[4]) # Consumer production efficiency
K = float(sys.argv[5]) # Resource carrying capacity

# Now define time -- integrate from 0 to 15, using 1000 points:
t = sc.linspace(0, 15,  1000)

x0 = 10
y0 = 5
z0 = sc.array([x0, y0]) # initials conditions: 10 prey and 5 predators per unit area

pops, infodict = integrate.odeint(dCR_dt, z0, t, full_output=True)

infodict['message']     # >>> 'Integration successful.'

prey, predators = pops.T # What's this for?
f1 = p.figure() #Open empty figure object
p.plot(t, prey, 'g-', label='Resource density') # Plot
p.plot(t, predators  , 'b-', label='Consumer density')
p.grid()
# add parameters to plot
paramstring = "r = " + str(r) + "\na = " + str(a) + "\nz = " + str(z) + "\ne = " + str(e) + "\nK = " + str(K)
p.annotate(paramstring,(14.1,2.1))
p.legend(loc='best')
p.xlabel('Time')
p.ylabel('Population')
p.title('Consumer-Resource population dynamics')
p.show()
f1.savefig('../Results/prey_and_predators_1.pdf')  # Save figure
