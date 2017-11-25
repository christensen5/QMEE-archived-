"""
This script runs example 1.1.9 (Computing derivaties) in the FEniCS book, but
translated to Firedrake. In order to do so, it re-solves the previous problem from membrane.py as well, but without plotting.
"""

from firedrake import *
import numpy
import matplotlib.pyplot as plt

#Initialise R, x0, y0 and sigma
R = 0.3  #radius of domain
theta = 0.2
x0 = 0.6*R*cos(theta)
y0 = 0.6*R*sin(theta)
sigma = 0.0025  # set >=50 for verification on a "flat" pressure

#======================================================================================================================
# SOLVE PREVIOUS PROBLEM
#======================================================================================================================
#Create the mesh and define the function space
mesh = Mesh("/home/alexander/Documents/QMEE/Firedrake Learning/circle.msh")  # Here we import a circular mesh made in gmsh
V = FunctionSpace(mesh, "Lagrange", 1)

#Define the variational problem
w = TrialFunction(V)
v = TestFunction(V)
f = Expression("4 * exp(-0.5 * pow((R * x[0] - x0)/sigma , 2) - 0.5 * pow((R * x[1] * y0)/sigma, 2))",
              R=R, x0=x0, y0=y0, sigma=sigma)  # define defaults for the parameters inside the new expression object (it can't see their existing values unless they're passed to it directly)
f = interpolate(f,V)
a = inner(nabla_grad(w), nabla_grad(v))*dx
L = f*v*dx

#Define boundary condition
bc = DirichletBC(V, Constant(0), (1))

#Compute solution
w = Function(V)
solve(a == L, w, bcs = [bc])


#======================================================================================================================
# NOW SOLVE THE GRADIENT PROBLEM
#======================================================================================================================
V_g = VectorFunctionSpace(mesh, "Lagrange", 1)
omega = TrialFunction(V_g)
v = TestFunction(V_g)

a = inner(omega, v)*dx
L = inner(grad(w), v)*dx

omega = Function(V_g)
solve(a == L, omega)

plot(omega, plot3d=True)
plt.show(block=True)
