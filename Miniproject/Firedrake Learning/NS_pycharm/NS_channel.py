from firedrake import *
import matplotlib.pyplot as plt
import numpy as np

# Initial stuff
T = 10.0  # final time
num_steps = 500  # number of time steps
dt = T / num_steps  # time step size

# Create mesh and define function spaces
mesh = UnitSquareMesh(16, 16)
V = VectorFunctionSpace(mesh, 'P', 2)  # P = piecewise, 2 = quadratic
Q = FunctionSpace(mesh, 'P', 1)  # piecewise linear
# together these make the Taylor-Hood element
Z = V * Q

# Define trial and test functions
# up = TrialFunctions(Z)
up = Function(Z)
u, p = split(up)
v, q = TestFunctions(Z)


# Define functions for solutions at previous and current time steps
up_n = Function(Z)
up_ = Function(Z)
u_n, p_n = split(up_n)  # previous time step
u_, p_ = split(up_)  # current time step

# Define boundary conditions
bcu = DirichletBC(Z.sub(0), as_vector([0.0, 0.0]), (3, 4))  # no slip on walls
bcp = [DirichletBC(Z.sub(1), Constant(8.0), 1),  # inflow pressure of 8
       DirichletBC(Z.sub(1), Constant(0.0), 2)]  # outflow pressure of 0

# Define expressions used in variational forms
U = 0.5 * (u_n + u)
n = FacetNormal(mesh)
f = Constant((0, 0))
k = Constant(dt)
mu = Constant(1)
rho = Constant(1)


# Define strain-rate tensor
def epsilon(u):
    return sym(nabla_grad(u))


# Define stress tensor
def sigma(u, p):
    return 2 * mu * epsilon(u) - p * Identity(len(u))


# Define variational problem for step 1
F1 = rho * inner((u - u_n) / k, v) * dx \
     + rho * inner(dot(u_n, nabla_grad(u_n)), v) * dx \
     + inner(sigma(U, p_n), epsilon(v)) * dx \
     + inner(p_n * n, v) * ds - inner(mu * nabla_grad(U) * n, v) * ds \
     - inner(f, v) * dx
a1 = lhs(F1)
L1 = rhs(F1)

# Define variational problem for step 2
a2 = dot(nabla_grad(p), nabla_grad(q)) * dx
L2 = dot(nabla_grad(p_n), nabla_grad(q)) * dx - (1 / k) * div(u_) * q * dx

# Define variational problem for step 3
a3 = dot(u, v) * dx
L3 = dot(u_, v) * dx - k * dot(nabla_grad(p_ - p_n), v) * dx

# Assemble matrices
A1 = assemble(a1, bcs=bcu)
A2 = assemble(a2, bcs=bcp)
A3 = assemble(a3)

# # Apply boundary conditions to matrices
# bcu.apply(A1)
# bcp[0].apply(A2)
# bcp[1].apply(A2)

# Time-stepping
t = 0
for n in range(num_steps):
    # Update current time
    t += dt

    # Step 1: Tentative velocity step
    b1 = assemble(L1, bcs=bcu)
    solve(A1, u_, b1)

    # Step 2: Pressure correction step
    b2 = assemble(L2, bcs=bcp)
    solve(A2, p_.vector(), b2)

    # Step 3: Velocity correction step
    b3 = assemble(L3)
    solve(A3, u_.vector(), b3)

    # Plot solution
    plot(u_)

    # # Compute error
    # u_e = Expression(('4*x[1]*(1.0 - x[1])', '0'), degree=2)
    # u_e = interpolate(u_e, V)
    # error = np.abs(u_e.vector().array() - u_.vector().array()).max()
    # print('t = %.2f: error = %.3g' % (t, error))
    # print('max u:', u_.vector().array().max())

    # Update previous solution
    u_n.assign(u_)
    p_n.assign(p_)

    # Hold plot
    # interactive()
