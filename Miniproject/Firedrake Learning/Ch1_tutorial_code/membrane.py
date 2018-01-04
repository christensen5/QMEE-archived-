"""
This script runs example 1.1.7 (warping of an elastic circular membrane with radius R) in the FEniCS book, but
translated to Firedrake.
"""

from firedrake import *
import numpy
import matplotlib.pyplot as plt

#Initialise R, x0, y0 and sigma
R = 0.3  #radius of domain
theta = 0.2
x0 = 0.6*R*cos(theta)
y0 = 0.6*R*sin(theta)
sigma = 50  # set >=50 for verification on a "flat" pressure

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

#Plot the solution and mesh
plot(w, plot3d=True)
plot(mesh)
plt.show(block=True)
# Save solution to file
File("/home/alexander/Documents/QMEE/Firedrake Learning/Ch1_tutorial_saves/1.1.6_membrane.pvd").write(w)

#Verification
w_exact = Expression("1 - x[0]*x[0] - x[1]*x[1]")
w_e = interpolate(w_exact, V)
dev = numpy.abs(w_e.vector().array() - w.vector().array()).max()
plot(w_e, plot3d=True)
print("sigma = ", sigma, "\nmax deviation = ", dev)