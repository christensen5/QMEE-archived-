""" A discrete-time Lotka-Volterra Model simulated using scipy, including random Gaussian prey fluctuations."""
__author__ = 'Alexander Christensen'
__version__ = '0.0.1'

import sys
import scipy.stats
import scipy as sc
import scipy.integrate as integrate
import pylab as p  # Contains matplotlib for plotting
from matplotlib.offsetbox import AnchoredText


def CR_t1(pops):
    """ Returns the growth rate of predator and prey populations at time step t1 given the previous time step t0."""

    # If first iteration, pops is only 1-dimensional
    # if pops.ndim == 1:
    #     R0 = pops[0] # find t0 Resource population
    #     C0 = pops[1] # find t0 Consumer population
    # else:
    R0 = pops[-1, 0] # find t0 Resource population
    C0 = pops[-1, 1] # find t0 Consumer population

    # Apply discrete LV equations to compute t1 Resource and Consumer populations
    R1 = R0 * (1 + (r + gscale * scipy.stats.norm.rvs()) * (1 - (R0 / K)) - a*C0)
    C1 = C0 * (1 - z + e * a * R0)

    return sc.array([R1, C1])

# Store parameters from user input:
r = float(sys.argv[1]) # Resource growth rate
a = float(sys.argv[2]) # Consumer search rate (determines consumption rate)
z = float(sys.argv[3]) # Consumer mortality rate
e = float(sys.argv[4]) # Consumer production efficiency
K = float(sys.argv[5]) # Resource carrying capacity
gscale = 0.125 # Scaling factor for random Gaussian component of growth rate.

# Now define time -- integrate from 0 to 30, using 1000 points:
t = sc.arange(0, 100) # CHANGING START POINT WILL NOT CHANGE START TIME OF SIMULATION BELOW!!

x0 = 10
y0 = 5
z0 = sc.array([x0, y0])  # initials conditions: 10 prey and 5 predators per unit area

pops = sc.array([[x0, y0]])
for timestep in t[:-1]: # All but the last timestep or we'll have one more population step than time steps
    pops = sc.append(pops, [CR_t1(pops)], axis=0)

prey = pops[:,0]
predators = pops[:,1]
f1 = p.figure() #Open empty figure object
p.plot(t, prey, 'g-', label='Resource density') # Plot
p.plot(t, predators, 'b-', label='Consumer density')
p.grid()
# add parameters to plot
paramstring = "r = " + str(r) + "\na = " + str(a) + "\nz = " + str(z) + "\ne = " + str(e) + "\nK = " + str(K) + "\neps ~ N(0," + str(gscale) + ")"
p.annotate(paramstring, (t.size - 10, 5.1))
p.legend(loc='best')
p.xlabel('Timesteps')
p.ylabel('Population')
p.title('Consumer-Resource discrete population dynamics (w/ random Prey fluctuation)')
p.show()
f1.savefig('../Results/prey_and_predators_5.pdf')  # Save figure
