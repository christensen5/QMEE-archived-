# In this script we will solve the 2D Navier-Stokes equation on a rectangular channel, using a
# splitting method as set out in the FEniCS example at https://fenicsproject.org/pub/tutorial/html/._ftut1009.html

# imports
from firedrake import *
import matplotlib.pyplot as plt
import numpy as np

# Initial settings
N = 16
T = 10.0  # final time
num_steps = 500  # time steps
dt = T/num_steps
mu_val = 1
rho_val = 1
outfile = File("/home/alexander/Documents/QMEE/Firedrake Learning/NS_tutorial_saves/NS_channel.pvd")

# Create mesh and function spaces for Step1
mesh = RectangleMesh(2*N, N, 2, 1)
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)
Z = V * Q

# Define test and trial functions
u, p = TrialFunctions(Z)
#u, p = split(up)
v, q = TestFunctions(Z)
#v, q = split(vq)


# Define functions for solutions at different time steps
x, y = SpatialCoordinate(mesh)
up_now = Function(Z, name="up_now")
up_next = Function(Z, name="up_next")
up_star = Function(Z, name="up_star")
u_now = Function(V, name="u_now").interpolate(as_vector((-x, 0*y)))
u_next = Function(V, name="u_next")
u_star = Function(V, name="u_star")
p_now = Function(Q, name="p_now")
p_next = Function(Q, name="p_next")


# Define expressions used in the variational forms
n = FacetNormal(mesh)
f = Constant((0.0, 0.0))
mu = Constant(mu_val)
rho = Constant(rho_val)
u_mid = 0.5*(u_now + u_next)
def epsilon(u):
    return sym(nabla_grad(u))
def sigma(u, p):
    return 2*mu*epsilon(u) - p*Identity(len(u))


# Define boundary conditions
bcu = DirichletBC(Z.sub(0), as_vector([0.0, 0.0]), (3, 4))  # no slip on walls
bcp = [DirichletBC(Z.sub(1), Constant(8.0), 1),  # inflow pressure of 8
       DirichletBC(Z.sub(1), Constant(0.0), 2)]  # outflow pressure of 0


# Define variational forms
F1 = inner(rho*(u - u_now)/dt, v) * dx \
    + inner(dot(rho*u_now, grad(u_now)), v) * dx \
    + inner(sigma(u_mid, p_now), epsilon(v)) * dx \
    + inner(p_now * n, v) * ds \
    - inner(mu * dot(grad(u_mid), n), v) * ds \
    - inner(f, v) * dx
a1, L1 = system(F1)

a2 = inner(grad(p), grad(q)) * dx
L2 = inner(grad(p_now), grad(q)) * dx \
    - (1/dt) * inner(div(u_star), q) * dx

a3 = inner(u, v) * dx
L3 = inner(u_star, v) * dx \
     - dt * inner(grad(p_next - p_now), v) * dx


# Time loop
t = 0.0
steps = 0
while t < T:
    # solve variational problems, extract components & increment t
    prob1 = LinearVariationalProblem(a1, L1, up_star, bcs=bcu)  # solve for u, will store as u_star
    solve1 = LinearVariationalSolver(prob1)
    solve1.solve()
    u_star, p_star = split(up_star)

    prob2 = LinearVariationalProblem(a2, L2, up_next, bcs=bcp)  # solve for p, will store as p_next
    solve2 = LinearVariationalSolver(prob2)
    solve2.solve()
    u_next, p_next = split(up_next)

    prob3 = LinearVariationalProblem(a3, L3, up_next)  # solve for u, will store as u_next
    solve3 = LinearVariationalSolver(prob3)
    solve3.solve()
    u_next, p_next = up_next.split()

    t += dt
    steps += 1

    # update solutions
    u_mid = (0.5*(u_now + u_next))  # Does this work??
    u_now.assign(u_next)
    p_now.assign(p_next)


    # write timestep and save solution every 50 steps
    if steps%50 == 0:
        outfile.write(u_next)
        print("t=", t)