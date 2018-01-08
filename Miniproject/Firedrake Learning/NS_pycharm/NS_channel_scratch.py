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
mu_val = 1 #0.001
rho_val = 1
outfile_u = File("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/NS_tutorial_saves/NS_channel/firedrake/u.pvd")
outfile_p = File("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/NS_tutorial_saves/NS_channel/firedrake/p.pvd")

# Create mesh and function spaces for Step1
mesh = UnitSquareMesh(16,16)
#mesh = RectangleMesh(2*N, N, 2, 1)
V = VectorFunctionSpace(mesh, "P", 2)
Q = FunctionSpace(mesh, "P", 1)
Z = V * Q

# Define test and trial functions
u, p = TrialFunctions(Z)
v, q = TestFunctions(Z)


# Define functions for solutions at different time steps
x, y = SpatialCoordinate(mesh)
up_now = Function(Z, name="up_now")
up_next = Function(Z, name="up_next")
up_star = Function(Z, name="up_star")
up_temp = Function(Z, name="up_temp")  # to hold solutions at each intermediate step
u_now = Function(V, name="u_now")#.interpolate(as_vector([y * (1 - y), 0]))
u_next = Function(V, name="u_next")#.interpolate(as_vector([y * (1 - y), 0]))
u_star = Function(V, name="u_star")
u_temp = Function(V, name="u_temp")  # to hold solutions at each intermediate step
p_now = Function(Q, name="p_now")
p_init = Function(Q).assign(p_now)#.interpolate(32 * y * (1 - y))  # save initial pressure
p_next = Function(Q, name="p_next")
p_temp = Function(Q, name="p_temp")  # to hold solutions at each intermediate step

# Save initial velocity and pressure states
outfile_u.write(u_next)
outfile_p.write(p_next)

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
bcu = DirichletBC(Z.sub(0), Constant(as_vector([0.0, 0.0])), (3, 4))  # no slip on walls
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


# Define linear problems
prob1 = LinearVariationalProblem(a1, L1, up_temp, bcs=bcu)  # solve for u, will store as up_star
solve1 = LinearVariationalSolver(prob1)#, solver_params={'ksp_type': 'bicgstab', 'pc_type': 'lu', 'pc_factor_mat_solver_type': 'umfpack',
                                                      # 'ksp_converged_reason': True,
                                                      # 'ksp_monitor_true_residual': True,
                                                      # 'ksp_view': True})

prob2 = LinearVariationalProblem(a2, L2, up_temp, bcs=bcp)  # solve for p, will store as up_next
solve2 = LinearVariationalSolver(prob2)#, solver_params={'ksp_type': 'bicgstab', 'pc_type': 'hypre_amg'})

prob3 = LinearVariationalProblem(a3, L3, up_temp)  # solve for u, will store as up_next
solve3 = LinearVariationalSolver(prob3)#, solver_params={'ksp_type': 'cg', 'pc_type': 'sor'})


# Time loop
t = 0.0
steps = 0
while t < T:
    # solve variational problems, extract components & increment t
    solve1.solve()
    u_temp, p_temp = up_temp.split()
    u_star.assign(u_temp)


    solve2.solve()
    u_temp, p_temp = up_temp.split()
    p_next.assign(p_temp)


    solve3.solve()
    u_temp, p_temp = up_temp.split()
    u_next.assign(u_temp)

    t += dt
    steps += 1

    # write timestep and save solution every 50 steps
    if steps%1 == 0:
        outfile_u.write(u_next)
        outfile_p.write(p_next)
        print("t=", t)

    # update solutions
    u_now.assign(u_next)
    p_now.assign(p_next)



