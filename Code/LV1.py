""" The typical Lotka-Volterra Model simulated using scipy """

import scipy as sc 
import scipy.integrate as integrate

import pylab as p #Contains matplotlib for plotting

# import matplotlip.pylab as p #Some people might need to do this

def dx_dt(z, t=0):
    """ Returns the growth rate of predator and prey populations at any 
    given time step """
    
    x = z[0]
    y = z[1]
    
    dxdt = r*x - a*x*y 
    dydt = -m*y + e*a*x*y
    
    return sc.array([dxdt, dydt])

# Define parameters:
r = 1. # Prey growth rate
a = 0.1 # Predator search rate (determines consumption rate) 
m = 1.5 # Predator mortality rate
e = 0.75 # Predator production efficiency

# Now define time -- integrate from 0 to 15, using 1000 points:
t = sc.linspace(0, 15,  1000)

x0 = 10
y0 = 5 
z0 = sc.array([x0, y0]) # initials conditions: 10 prey and 5 predators per unit area

z, infodict = integrate.odeint(dx_dt, z0, t, full_output=True)

infodict['message']                     # >>> 'Integration successful.'

prey, predators = z.T # What's this for?
f1 = p.figure() #Open empty figure object
p.plot(t, prey, 'g-', label='Prey density') # Plot
p.plot(t, predators  , 'b-', label='Predator density')
p.grid()
p.legend(loc='best')
p.xlabel('Time')
p.ylabel('Population')
p.title('Predator-prey population dynamics')
p.show()
f1.savefig('prey_and_predators_1.pdf') #Save figure
