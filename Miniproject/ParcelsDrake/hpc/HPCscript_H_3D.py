from incflow import *
from incflow.inc_navier_stokes_3D import IncNavierStokes3D
from firedrake import *
import time, os

start = time.time()
walltime = 300  # in seconds
savetime = 30

mesh = Mesh("NavierStokes/meshes/H_3D_mixed.msh")
rho = 1
nu = 0.01
dt = 0.5

INS = IncNavierStokes3D(mesh, nu, rho, dt, verbose=True)
W = INS.get_mixed_fs()
u_bcs = [DirichletBC(W.sub(0), Constant((0.0, 0.0, 0.0)), 70),
         DirichletBC(W.sub(0), Constant((0.0, 1.0, 0.0)), 68),  # inflow velocity of (0, 1, 0)
         DirichletBC(W.sub(0), Constant((0.0, 1.0, 0.0)), 69)]  # outflow velocity of (0, 1, 0)
p_bcs = [DirichletBC(W.sub(1), Constant(0.0), 69)]  # outflow pressure of 0
INS.set_bcs(u_bcs, p_bcs)
INS.setup_solver()

step = 0
t = 0

folderstr = str(t) + "Ti_" + str(t_end) + "Tf_" + str(dt) + "dt_" + str(rho*nu) + "mu"
outdir = os.getcwd() + "/data/" + folderstr
outfile_u = File(outdir + "/u.pvd")
chk_out = checkpointing.HDF5File(outdir + "/u_field_0.h5", file_mode='w')

while time.time()-start < walltime - savetime:
    step += 1
    t += dt
    u_sol, p_sol = INS.step()
    chk_out.new_file(name="u_field_" + str(step))
    chk_out.store(u_sol)

outfile_u.write(u_sol)
print("End simulation time = " + str(t))
