from incflow import *
from incflow.inc_navier_stokes_3D import IncNavierStokes3D
from firedrake import *
from tqdm import tqdm

mesh = Mesh("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/meshes/H_3D_mixed.msh")
rho = 1
nu = 0.01
dt = 0.01

INS = IncNavierStokes3D(mesh, nu, rho, dt, verbose=True)
W = INS.get_mixed_fs()
u_bcs = [DirichletBC(W.sub(0), Constant((0.0, 0.0, 0.0)), 70)]
p_bcs = [DirichletBC(W.sub(1), Constant(8.0), 68),  # inflow pressure of 8
           DirichletBC(W.sub(1), Constant(0.0), 69)]  # outflow pressure of 0
INS.set_bcs(u_bcs, p_bcs)
INS.setup_solver()

step = 0
t = 0.0
t_end = 10
num_steps = int((t_end - t)/INS.dt)
folderstr = str(t_end) + "T_" + str(dt) + "dt_" + str(rho*nu) + "mu"
outfile_u = File("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/NS_tutorial_saves/NS_H_3D/Julian/" + folderstr + "/u.pvd")
outfile_p = File("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/NS_tutorial_saves/NS_H_3D/Julian/" + folderstr + "/p.pvd")

for steps in tqdm(range(num_steps)):
    t += INS.dt

    u_sol, p_sol = INS.step()

    if steps % 10 == 0:
        outfile_u.write(u_sol)
        outfile_p.write(p_sol)