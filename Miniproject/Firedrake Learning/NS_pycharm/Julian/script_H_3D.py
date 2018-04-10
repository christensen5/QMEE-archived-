from incflow import *
from incflow.inc_navier_stokes_3D import IncNavierStokes3D
from firedrake import *
from tqdm import tqdm

mesh = Mesh("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/meshes/H_3D_mixed.msh")
rho = 1
nu = 0.01
dt = 0.5

INS = IncNavierStokes3D(mesh, nu, rho, dt, verbose=True)
W = INS.get_mixed_fs()
# u_bcs = [DirichletBC(W.sub(0), Constant((0.0, 0.0, 0.0)), 70)]
# p_bcs = [DirichletBC(W.sub(1), Constant(8.0), 68),  # inflow pressure of 8
#            DirichletBC(W.sub(1), Constant(0.0), 69)]  # outflow pressure of 0
u_bcs = [DirichletBC(W.sub(0), Constant((0.0, 0.0, 0.0)), 70),
         DirichletBC(W.sub(0), Constant((0.0, 1.0, 0.0)), 68),  # inflow velocity of (0, 1, 0)
         DirichletBC(W.sub(0), Constant((0.0, 1.0, 0.0)), 69)]  # outflow velocity of (0, 1, 0)
p_bcs = [DirichletBC(W.sub(1), Constant(0.0), 69)]  # outflow pressure of 0
INS.set_bcs(u_bcs, p_bcs)
#up_init_path = "/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/inputs/Hmixed/data1/0.5dt_0.01mu/h5/u_func_01437.h5"
INS.setup_solver()#up_init=up_init_path)

step = 0
t = 0
t_end = 600
num_steps = int((t_end - t)/INS.dt)
folderstr = str(t) + "Ti_" + str(t_end) + "Tf_" + str(dt) + "dt_" + str(rho*nu) + "mu" + "_1cpu"
outfile_u = File("/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/inputs/" + folderstr + "/u.pvd")
#outfile_p = File("/media/alexander/DATA/Ubuntu/Miniproject/Firedrake Learning/outputs/NS_H_3D/Julian/" + folderstr + "/p.pvd")

for steps in tqdm(range(num_steps)):
    t += INS.dt

    u_sol, p_sol = INS.step()

    if steps % 10 == 0:
        outfile_u.write(u_sol)
        #outfile_p.write(p_sol)



# Checkpoint final state.
chk_out = checkpointing.HDF5File("/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/inputs/" + folderstr + "/final.h5", file_mode='w')
chk_out.write(u_sol, "/velocity")
chk_out.write(p_sol, "/pressure")
chk_out.write(INS.up, "/up")
chk_out.close()