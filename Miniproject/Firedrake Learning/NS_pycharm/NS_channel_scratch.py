from firedrake import *

# Hard constants
mu_val = 1
rho_val = 1
T = 10.0           # final time
num_steps = 500    # number of time steps
dt = T / num_steps # time step size

# Mesh, functionspaces and functions
mesh = UnitSquareMesh(16,16)
V = VectorFunctionSpace(mesh, "P", 2)
Q = FunctionSpace(mesh, "P", 1)
u = TrialFunction(V)
v = TestFunction(V)
p = TrialFunction(V)
q = TestFunction(Q)
u_now = Function(V)
u_next = Function(V)
u_star = Function(V)
p_now = Function(Q)
p_next = Function(Q)

# Expressions for the variational forms
n = FacetNormal(mesh)
f = Constant((0.0, 0.0))
k = Constant(dt)
mu = Constant(mu_val)
rho = Constant(mu_val)
u_mid = 0.5*(u_now + u_next)
def sigma(u, p):
    return 2*mu*sym(nabla_grad(u)) - p*Identity(len(u))

# Define boundary conditions
bcu = DirichletBC(V, Constant([0.0, 0.0]), (3, 4))  # no slip on walls
bcp = [DirichletBC(Q, Constant(8.0), 1),  # inflow pressure of 8
       DirichletBC(Q, Constant(0.0), 2)]  # outflow pressure of 0

# Define variational forms
F1 = inner(rho*(u - u_now)/k, v) * dx \
    + inner(dot(rho*u_now, nabla_grad(u_now)), v) * dx \
    + inner(sigma(u_mid, p_now), sym(nabla_grad(v))) * dx \
    + inner(p_now * n, v) * ds \
    - inner(mu * dot(nabla_grad(u_mid), n), v) * ds \
    - inner(f, v) * dx
a1, L1 = system(F1)

a2 = inner(nabla_grad(p), nabla_grad(q)) * dx
L2 = inner(nabla_grad(p_now), nabla_grad(q)) * dx \
    - (1/k) * inner(div(u_star), q) * dx

a3 = inner(u, v) * dx
L3 = inner(u_star, v) * dx \
     - k * inner(nabla_grad(p_next - p_now), v) * dx


# Define linear problems
prob1 = LinearVariationalProblem(a1, L1, u_star, bcs=bcu)
prob2 = LinearVariationalProblem(a2, L2, p_next, bcs=bcp)
prob3 = LinearVariationalProblem(a3, L3, u_next)

# Time loop
t = 0.0
steps = 0
while t < 10.0:
    solve1 = LinearVariationalSolver(prob1)
    solve1.solve()

    solve2 = LinearVariationalSolver(prob2)
    solve2.solve()

    solve3 = LinearVariationalSolver(prob3)
    solve3.solve()

    t += dt
    steps += 1

    # update solutions
    u_now.assign(u_next)
    p_now.assign(p_next)

