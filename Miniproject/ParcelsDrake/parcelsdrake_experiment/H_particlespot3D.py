# Imports
from datetime import timedelta

import numpy as np
from firedrake import *
from parcels import *
from tqdm import tqdm

from incflow.inc_navier_stokes_3D import IncNavierStokes3D
from pdist import uniform_H_Dist

def extract_field_3D(u, grid_start, grid_end, grid_incr, depth):
    lon = np.arange(grid_start, grid_end, grid_incr)
    lat = np.arange(grid_start, grid_end, grid_incr)
    dep = np.arange(0, depth, grid_incr)
    u_field = np.empty([len(lon), len(lat), len(dep), 3])
    for x in range(len(lon)):
        for y in range(len(lat)):
            for z in range(len(dep)):
                try:
                    u_field[x, y, z] = u.at([lon[x], lat[y], dep[z]])
                except PointNotInDomainError:
                    u_field[x, y, z] = [0, 0, 0]
    return u_field


# Navier-Stokes firedrake setup
mesh_NS = Mesh("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/meshes/H_3D_low.msh")
rho = 1
nu = 0.01
dt_NS = 0.1
INS = IncNavierStokes3D(mesh_NS, nu, rho, dt_NS)
W = INS.get_mixed_fs()
u_bcs = [DirichletBC(W.sub(0), Constant((0.0, 0.0, 0.0)), 70),
         DirichletBC(W.sub(0), Constant((0.0, 1.0, 0.0)), 68),  # inflow velocity of (0, 1, 0)
         DirichletBC(W.sub(0), Constant((0.0, 1.0, 0.0)), 69)]  # outflow velocity of (0, 1, 0)
p_bcs = [DirichletBC(W.sub(1), Constant(0.0), 69)]  # outflow pressure of 0
INS.set_bcs(u_bcs, p_bcs)
#up_init_path = "/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/inputs/0Ti_30Tf_0.1dt_0.01mu_velBCs.h5"
INS.setup_solver()#up_init=up_init_path)
u_sol, p_sol = INS.up.split()


# Lagrangian parcels setup
num_particles = 500
init_particles = 100
max_replace = 50
repeat_particles = 100
grid_res = 100
grid_start = 0
grid_end = 13.6
depth = 2.8
grid_incr = abs(grid_start - grid_end)/grid_res
lon = np.arange(grid_start, grid_end, grid_incr)
lat = np.arange(grid_start, grid_end, grid_incr)
depth_array = np.array([[np.arange(0, depth, grid_incr) for x in lon] for y in lat])
grid = RectilinearSGrid(lon=lon, lat=lat, depth=depth_array, time=np.zeros(1), mesh='flat')

Ufield_init = Field(name='U', data=np.zeros((len(lon), len(lat), int(np.floor(depth/grid_incr)+1))), grid=grid, allow_time_extrapolation=True)
Vfield_init = Field(name='V', data=np.zeros((len(lon), len(lat), int(np.floor(depth/grid_incr)+1))), grid=grid, allow_time_extrapolation=True)
#Wfield_init = Field(name='W', data=np.zeros((len(lon), len(lat), depth/grid_incr)), grid=grid, allow_time_extrapolation=True)
fieldset = FieldSet(U=Ufield_init, V=Vfield_init)

class FiredrakeParticle3D(JITParticle):
    u = Variable('u', dtype=np.float32, initial=0.)
    v = Variable('v', dtype=np.float32, initial=0.)
    w = Variable('w', dtype=np.float32, initial=0.)
    # vort_i = Variable('vort_i', dtype=np.float32, initial=0.)
    # vort_j = Variable('vort_j', dtype=np.float32, initial=0.)
    # vort_k = Variable('vort_k', dtype=np.float32, initial=0.)
    # dir_x = Variable('dir_x', dtype=np.float32, initial=0.)
    # dir_y = Variable('dir_y', dtype=np.float32, initial=0.)
    # dir_z = Variable('dir_z', dtype=np.float32, initial=0.)

repeatdt = 0.5
pset = ParticleSet(fieldset=fieldset, pclass=FiredrakeParticle3D, lon=np.repeat(2.8, repeat_particles), lat=np.linspace(4.8, 8.8, repeat_particles), depth=np.repeat(1.4, repeat_particles))#, repeatdt=repeatdt)

# Custom Parcels kernels
def DeleteParticle(particle, fieldset, time, dt):  # delete particles who run out of bounds.
    particle.delete()


# Time loop setup
t = 0
t_add_particles = 60
t_end = 120
parcels_interval = 1
num_steps = int((t_end - t)/INS.dt)
folderstr = str(t) + "Ti_" + str(t_end) + "Tf_" + str(dt_NS) + "dt_" + str(rho*nu) + "mu"
outfile_u = File("/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/outputs/H_3D/particlespot/" + folderstr + "/u.pvd")
outfile_particles = pset.ParticleFile(name="/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/outputs/H_3D/particlespot/" + folderstr + "/3D_test", outputdt=timedelta(seconds=dt_NS))#, type='indexed')
for steps in tqdm(range(num_steps)):
    t += dt_NS

    # Update firedrake solution
    u_sol, p_sol = INS.step()
    #u_vort = curl(u_sol)

    # Checkpoint final pre-particle state
    if t == t_add_particles:
        chk_out = checkpointing.HDF5File(
            "/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/outputs/H_3D/fields/" + str(t) + "s_" + str(
                rho * nu) + "mu_" + "Hlow_" + "preparticle.h5", file_mode='w')
        chk_out.write(u_sol, "/velocity")
        chk_out.write(p_sol, "/pressure")
        chk_out.write(INS.up, "/up")
        chk_out.close()

    # Update parcels field for plotting (slow so dont do often!).
    if t >= t_add_particles and steps % 10 == 0:
        field_data = extract_field_3D(u_sol, grid_start=grid_start, grid_end=grid_end, grid_incr=grid_incr, depth=depth)
        # plotfieldset = FieldSet.from_data(data={'U': field_data[:, :, 0], 'V': field_data[:, :, 1]},
        #                                   dimensions={'lon': np.arange(grid_start, grid_end - grid_incr, grid_incr),
        #                                               'lat': np.arange(grid_start, grid_end, grid_incr),
        #                                               'time': np.array([t], float)},
        #                                   transpose=False, mesh='flat')

    if (t >= t_add_particles):  # and (steps%parcels_interval==0):
        # Extract velocity field at particle positions
        for particle in pset:
            try:
                particle.u, particle.v, particle.w = u_sol.at([particle.lon, particle.lat, particle.depth])
            except PointNotInDomainError:
                tqdm.write("\nRemoved out-of-domain particle %d at point (%f, %f, %f) during velocity extraction." % (
                particle.id, particle.lon, particle.lat, particle.depth))
                particle.delete()
                continue

        # Update particle positions
        pset.execute(AdvectionEE_firedrake_3D,  # the kernel (which defines how particles move)
                     runtime=timedelta(seconds=dt_NS * parcels_interval),  # the total length of the run
                     dt=timedelta(seconds=dt_NS),  # the timestep of the kernel
                     recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle},
                     output_file=outfile_particles)

        # Save image to disk
        outfile_u.write(u_sol)


# Checkpoint final state.
chk_out = checkpointing.HDF5File(
    "/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/outputs/H_3D/fields/" + str(t_end) + "s_" + str(
        rho * nu) + "mu_" + "Hlow_" + "final.h5", file_mode='w')
chk_out.write(u_sol, "/velocity")
chk_out.write(p_sol, "/pressure")
chk_out.write(INS.up, "/up")
chk_out.close()

print("Ended")