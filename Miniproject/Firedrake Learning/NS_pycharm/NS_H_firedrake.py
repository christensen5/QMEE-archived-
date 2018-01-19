import sys
from firedrake import *
import matplotlib.pyplot as plt
from tqdm import tqdm  # progress bar


def NS_H_firedrake(viscosity=0.001, T=0.5, num_steps=5000, save_interval=50, u_init=False, p_init=False, bootstrap=False, solver_params1={}, solver_params2={}, solver_params3={}):

    # Hard constants
    density = 1
    dt = T / num_steps  # time step size
    if save_interval == "firstlast":
        save_interval = num_steps-2

    # Mesh, function spaces and functions
    #mesh = Mesh("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/meshes/Trot.msh")
    mesh = Mesh("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/meshes/H.msh")
    V = VectorFunctionSpace(mesh, "P", 2)
    V_diff = VectorFunctionSpace(mesh, 'DG', 1)
    Q = FunctionSpace(mesh, "P", 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    p = TrialFunction(Q)
    q = TestFunction(Q)
    if u_init:
        u_now = Function(V)
        chk_in = checkpointing.HDF5File("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/NS_tutorial_saves/NS_H/firedrake/dump.h5", file_mode='r')
        chk_in.read(u_now, "/velocity")
        chk_in.close()
    else:
        u_now = Function(V)
    u_next = Function(V)
    u_star = Function(V)
    if p_init:
        p_now = Function(Q)
        chk_in = checkpointing.HDF5File("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/NS_tutorial_saves/NS_H/firedrake/dump.h5", file_mode='r')
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
    ######################## Trot BCs ##################################
    # bcu = DirichletBC(V, Constant((0.0, 0.0)), 20)  # no slip on walls
    # bcp = [DirichletBC(Q, Constant(8.0), 21),  # inflow pressure of 8
    #        DirichletBC(Q, Constant(0.0), 22)]  # outflow pressure of 0
    ########################   H BCs  ##################################

    u_diff = Function(V_diff).interpolate(Dx(u_star, 1))
    bcu = [DirichletBC(V, Constant((0.0, 0.0)), 15),  # no slip on walls
           DirichletBC(V, Constant((0.0, 1.0)), 16)]  # inflow velocity of (0,1)
    bcp = DirichletBC(Q, u_diff[1], 17)  # outflow pressure of du_2/dx_2 (but index like a programmer not a mathmo)
                                                # (and use u_star since this is applied in the second linear problem, which computes p using u_star)

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
    outfile_u = File("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/NS_tutorial_saves/NS_H/firedrake/u.pvd")
    outfile_p = File("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/NS_tutorial_saves/NS_H/firedrake/p.pvd")
    outfile_u.write(u_save)
    outfile_p.write(p_save)

    # Time loop
    t = 0.0
    for steps in tqdm(range(num_steps)):
        solve1 = LinearVariationalSolver(prob1, solver_parameters=solver_params1)
        solve1.solve()

        solve2 = LinearVariationalSolver(prob2, solver_parameters=solver_params2)
        solve2.solve()

        solve3 = LinearVariationalSolver(prob3, solver_parameters=solver_params3)#, 'pc_type': 'lu', 'pc_factor_mat_solver_package': 'mumps'})
        solve3.solve()

        t += dt

        if steps%save_interval==1 or steps==num_steps:
            u_save.assign(u_next)
            p_save.assign(p_next)
            outfile_u.write(u_save)
            outfile_p.write(p_save)

        # update solutions
        u_now.assign(u_next)
        p_now.assign(p_next)

    # If running in bootstrap mode, return the final solutions for u and p to be used as initial conditions for the next
    # run.
    if bootstrap:
        chk_out = checkpointing.HDF5File("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/NS_tutorial_saves/NS_H/firedrake/dump.h5", file_mode='w')
        chk_out.write(u_now, "/velocity")
        chk_out.write(p_now, "/pressure")
        chk_out.close()

def bootstrapper():
    solver_params1 = {'ksp_type': 'bcgs', 'pc_type': 'hypre'}
    solver_params2 = {'ksp_type': 'bicg', 'pc_type': 'hypre'}
    solver_params3 = {'ksp_type': 'cg', 'pc_type': 'sor'}
    NS_H_firedrake(1, 5, 5000, "firstlast", False, False, True, solver_params1, solver_params2, solver_params3)
    tqdm.write("mu=1 bootstrapping complete!")
    NS_H_firedrake(0.1, 20, 20000, "firstlast", True, True, True, solver_params1, solver_params2, solver_params3)
    tqdm.write("mu=0.1 bootstrapping complete!")
    NS_H_firedrake(0.01, 5, 10000, 500, True, True, True, solver_params1, solver_params2, solver_params3)
    tqdm.write("mu=0.01 bootstrapping complete!")
    NS_H_firedrake(0.001, 5, 10000, 100, True, True, False, solver_params1, solver_params2, solver_params3)


if __name__ == "__main__":
    solver_params1 = {'ksp_type': 'bcgs', 'pc_type': 'hypre'}
    solver_params2 = {'ksp_type': 'bicg', 'pc_type': 'hypre'}
    solver_params3 = {'ksp_type': 'cg', 'pc_type': 'sor'}
    NS_H_firedrake(1, 5, 5000, 100, False, False, True, solver_params1, solver_params2, solver_params3)
    #bootstrapper()
    input("Press Enter to end.")