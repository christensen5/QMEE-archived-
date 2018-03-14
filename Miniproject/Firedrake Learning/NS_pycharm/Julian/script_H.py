import os, shutil
from incflow import *
from incflow.inc_navier_stokes import IncNavierStokes
from firedrake import *
from tqdm import tqdm

mesh = Mesh("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/meshes/H.msh")
rho = 1
nu = 0.001
dt = 0.1

INS = IncNavierStokes(mesh, nu, rho, dt)
W = INS.get_mixed_fs()
#u_bcs = [DirichletBC(W.sub(0), Constant((0.0, 0.0)), 15)]
#p_bcs = [DirichletBC(W.sub(1), Constant(8.0), 16),  # inflow pressure of 8
#           DirichletBC(W.sub(1), Constant(0.0), 17)]  # outflow pressure of 0
u_bcs = [DirichletBC(W.sub(0), Constant((0.0, 0.0)), 15),
         DirichletBC(W.sub(0), Constant((0.0, 1.0)), 16),  # inflow velocity of (0,1)
         DirichletBC(W.sub(0), Constant((0.0, 1.0)), 17)]  # outflow velocity of (0,1)
p_bcs = [DirichletBC(W.sub(1), Constant(0.0), 17)]  # outflow pressure of 0

INS.set_bcs(u_bcs, p_bcs)
INS.setup_solver()

step = 0
t = 0.0
t_end = 60
num_steps = int((t_end - t)/INS.dt)

simstr = str(t) + "Ti_" + str(t_end) + "Tf_" + str(dt) + "dt_" + str(rho*nu) + "mu"
folderstr = "/media/alexander/DATA/Ubuntu/Miniproject/Firedrake Learning/outputs/NS_H/Julian/" + simstr + "_velBC"
# if os.path.isdir(folderstr):  # Empty directory if already existing.
#     shutil.rmtree(folderstr)
#     os.mkdir(folderstr)
outfile_u = File(folderstr + "/u.pvd")
outfile_p = File(folderstr + "/p.pvd")

for steps in tqdm(range(num_steps)):
    t += INS.dt

    u_sol, p_sol = INS.step()

#    if steps % 10 == 0:
    outfile_u.write(u_sol)
    outfile_p.write(p_sol)


# Checkpoint final state.
chk_out = checkpointing.HDF5File(folderstr + "/60s_0.001mu_H.h5", file_mode='w')
chk_out.write(u_sol, "/velocity")
chk_out.write(p_sol, "/pressure")
chk_out.write(INS.up, "/up")
chk_out.write(t, "/t")
chk_out.close()

