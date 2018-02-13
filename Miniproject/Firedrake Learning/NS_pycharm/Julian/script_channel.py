from Julian import *
from Julian.inc_navier_stokes import IncNavierStokes
from firedrake import *
from tqdm import tqdm

mesh = RectangleMesh(32, 16, 2, 1)
rho = 1
nu = 0.1
dt = 0.001

outfile_u = File("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/NS_tutorial_saves/NS_channel/Julian/u.pvd")
outfile_p = File("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/NS_tutorial_saves/NS_channel/Julian/p.pvd")

INS = IncNavierStokes(mesh, nu, rho, dt)
W = INS.get_mixed_fs()
u_bcs = [DirichletBC(W.sub(0), Constant((0.0, 0.0)), (3, 4))]
p_bcs = [DirichletBC(W.sub(1), Constant(8.0), 1),  # inflow pressure of 8
           DirichletBC(W.sub(1), Constant(0.0), 2)]  # outflow pressure of 0
INS.set_bcs(u_bcs, p_bcs)
INS.setup_solver()

step = 0
t = 0.0
t_end = 1
num_steps = int((t_end - t)/INS.dt)


for steps in tqdm(range(num_steps)):
    t += INS.dt

    u_sol, p_sol = INS.step()

    if steps % 50 == 0:
        outfile_u.write(u_sol)
        outfile_p.write(p_sol)


