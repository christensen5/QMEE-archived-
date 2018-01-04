#[1]pg5

from firedrake import *

# Create a mesh and define function space
mesh = UnitSquareMesh(6,4)
V = FunctionSpace(mesh, "Lagrange", 1)

# Define boundary conditions
u0 = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]")

bc = DirichletBC(V, u0, (1, 2, 3, 4))

# Define the variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = inner(nabla_grad(u), nabla_grad(v))*dx
L = f*v*dx

# Compute solution
u = Function(V)
solve(a == L, u, bcs=[bc])

# Plot solution and mesh
plot(u)
plot(mesh)
