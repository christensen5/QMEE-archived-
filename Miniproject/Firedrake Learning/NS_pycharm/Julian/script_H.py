from incflow import *
from incflow.inc_navier_stokes import IncNavierStokes
from firedrake import *
from tqdm import tqdm

mesh = Mesh("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/meshes/H_low.msh")
rho = 1
nu = 0.01
dt = 0.05

INS = IncNavierStokes(mesh, nu, rho, dt)
W = INS.get_mixed_fs()
u_bcs = [DirichletBC(W.sub(0), Constant((0.0, 0.0)), 15)]
p_bcs = [DirichletBC(W.sub(1), Constant(8.0), 16),  # inflow pressure of 8
           DirichletBC(W.sub(1), Constant(0.0), 17)]  # outflow pressure of 0
INS.set_bcs(u_bcs, p_bcs)
INS.setup_solver()

step = 0
t = 0.0
t_end = 600
num_steps = int((t_end - t)/INS.dt)

folderstr = str(t_end) + "T_" + str(dt) + "dt_" + str(rho*nu) + "mu"
outfile_u = File("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/NS_tutorial_saves/NS_H/Julian/" + folderstr + "/u.pvd")
outfile_p = File("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/NS_tutorial_saves/NS_H/Julian/" + folderstr + "/p.pvd")

for steps in tqdm(range(num_steps)):
    t += INS.dt

    u_sol, p_sol = INS.step()

    if steps % 50 == 0:
        outfile_u.write(u_sol)
        outfile_p.write(p_sol)


