import sys
from firedrake import *
import matplotlib.pyplot as plt
import time
from tqdm import tqdm  # progress bar


def NS_cylinder(viscosity=0.001, T=0.5, num_steps=5000, save_interval=50, u_init=False, p_init=False, solver_params1={'ksp_type': 'bcgs', 'pc_type': 'hypre'}, solver_params2={'ksp_type': 'bcgs', 'pc_type': 'hypre'}, solver_params3={'ksp_type': 'cg', 'pc_type': 'sor'}):

    # Hard constants
    density = 1
    dt = T / num_steps  # time step size
    if save_interval == "firstlast":
        save_interval = num_steps-2

    # Mesh, function spaces and functions
    #mesh = Mesh("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/meshes/Trot.msh")
    mesh = Mesh("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/meshes/rectcylinder.msh")
    V = VectorFunctionSpace(mesh, "P", 2)
    V_xy = FunctionSpace(mesh, "P", 2)
    Q = FunctionSpace(mesh, "P", 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    p = TrialFunction(Q)
    q = TestFunction(Q)
    if u_init:
        u_now = Function(V)
        chk_in = checkpointing.HDF5File("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/NS_tutorial_saves/NS_cylinder/firedrake/init.h5", file_mode='r')
        chk_in.read(u_now, "/velocity")
        chk_in.close()
    else:
        u_now = Function(V)
    u_next = Function(V)
    u_star = Function(V)
    if p_init:
        p_now = Function(Q)
        chk_in = checkpointing.HDF5File("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/NS_tutorial_saves/NS_cylinder/firedrake/init.h5", file_mode='r')
        chk_in.read(p_now, "/pressure")
        chk_in.close()
    else:
        p_now = Function(Q)
    p_next = Function(Q)


    # Expressions for the variational forms
    n = FacetNormal(mesh)
    f = Constant((0.0, 0.0))
    k = Constant(dt)
    mu = Constant(viscosity)
    rho = Constant(density)
    u_mid = 0.5*(u_now + u)
    def sigma(u, p):
        return 2*mu*sym(nabla_grad(u)) - p*Identity(len(u))

    # Define boundary conditions
    bcu = DirichletBC(V, Constant((0.0, 0.0)), (1, 4))
    bcp = [DirichletBC(Q, Constant(8.0), 2),  # inflow pressure of 8
           DirichletBC(Q, Constant(0.0), 3)]  # outflow pressure of 0

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
    prob3 = LinearVariationalProblem(a3, L3, u_next, bcs=bcu)

    # Prep for saving solutions
    u_save = Function(V).assign(u_now)
    p_save = Function(Q).assign(p_now)
    outfile_u = File("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/NS_tutorial_saves/NS_cylinder/firedrake/u.pvd")
    outfile_p = File("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/NS_tutorial_saves/NS_cylinder/firedrake/p.pvd")
    outfile_u.write(u_save)
    outfile_p.write(p_save)

    # Time loop
    t = 0.0

    for steps in tqdm(range(num_steps)):
        solve1 = LinearVariationalSolver(prob1, solver_parameters=solver_params1)
        solve1.solve()

        solve2 = LinearVariationalSolver(prob2, solver_parameters=solver_params2)
        solve2.solve()

        solve3 = LinearVariationalSolver(prob3, solver_parameters=solver_params3)
        solve3.solve()

        t += dt

        # Save solutions and numpy velocity field arrays at appropriate timesteps
        if steps%save_interval==1 or steps==num_steps:
            u_save.assign(u_next)
            p_save.assign(p_next)
            outfile_u.write(u_save)
            outfile_p.write(p_save)


        # update solutions
        u_now.assign(u_next)
        p_now.assign(p_next)


if __name__ == "__main__":
    solver_params1 = {'ksp_type': 'bcgs', 'pc_type': 'hypre'}  # bgcs, hypre
    solver_params2 = {'ksp_type': 'bcgs', 'pc_type': 'hypre'}
    solver_params3 = {'ksp_type': 'bcgs', 'pc_type': 'hypre'}  # cg, sor
    NS_cylinder(0.01, 60, 60000, 60000, False, False, False, solver_params1, solver_params2, solver_params3)
    #input("Press Enter to end.")