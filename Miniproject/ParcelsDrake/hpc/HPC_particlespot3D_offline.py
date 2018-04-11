# Imports
from datetime import timedelta
import os

import numpy as np
from firedrake import *
from parcels import *

from incflow.inc_navier_stokes_3D import IncNavierStokes3D
from pdist import uniform_H_dist, uniform_H_bar_entry_dist

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
mesh_NS = Mesh("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/meshes/H_3D_mixed.msh")
rho = 1
nu = 0.01
dt_NS = 0.5
INS = IncNavierStokes3D(mesh_NS, nu, rho, dt_NS)
W = INS.get_mixed_fs()
up_sol = Function(W)
path_to_fields = "/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/inputs/Hmixed/data1/0.5dt_0.01mu/h5/"
chk_in = checkpointing.HDF5File(path_to_fields + "u_func_00001.h5", file_mode='r')
chk_in.read(up_sol, "/up")
chk_in.close()
u_sol, p_sol = up_sol.split()


# Lagrangian parcels setup
init_particles = 1000
grid_res = 100
grid_start = 0
grid_end = 13.6
depth = 2.8
grid_incr = abs(grid_start - grid_end)/grid_res
lon = np.arange(grid_start, grid_end, grid_incr)
lat = np.arange(grid_start, grid_end, grid_incr)
grid = RectilinearZGrid(lon=lon, lat=lat, depth=np.arange(0, depth, grid_incr), time=np.zeros(1), mesh='flat')

Ufield_init = Field(name='U', data=np.zeros((int(np.floor(depth/grid_incr) + 1), len(lat), len(lon))), grid=grid, allow_time_extrapolation=True)
Vfield_init = Field(name='V', data=np.zeros((int(np.floor(depth/grid_incr) + 1), len(lat), len(lon))), grid=grid, allow_time_extrapolation=True)
fieldset = FieldSet(U=Ufield_init, V=Vfield_init)

class FiredrakeParticle3D(JITParticle):
    u = Variable('u', dtype=np.float32, initial=0.)
    v = Variable('v', dtype=np.float32, initial=0.)
    w = Variable('w', dtype=np.float32, initial=0.)


# pfield_uniform = Field(name='pfield_uniform', data=uniform_H_dist(grid), transpose=False, grid=grid)
# pset = ParticleSet.from_field(fieldset=fieldset,
#                               pclass=FiredrakeParticle3D,
#                               start_field=pfield_uniform,
#                               depth=np.random.rand(init_particles) * depth,
#                               size=init_particles)
pfield_Hbar_entry = Field(name='pfield_Hbar_entry', data=uniform_H_bar_entry_dist(grid), transpose=False, grid=grid)
pset = ParticleSet.from_field(fieldset=fieldset,
                              pclass=FiredrakeParticle3D,
                              start_field=pfield_Hbar_entry,
                              depth=np.linspace(0, depth, init_particles),
                              size=init_particles)

# Custom Parcels kernels
def DeleteParticle(particle, fieldset, time, dt):  # delete particles who run out of bounds.
    particle.delete()


# Time loop setup
t = 0
t_add_particles = 60
t_end = 65
parcels_interval = 1
num_steps = int((t_end - t)/dt_NS)

outpath = "/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/outputs/H_3D/particlespot/offline/"
folderstr = str(t) + "Ti_" + str(t_end) + "Tf_" + str(dt_NS) + "dt_" + str(rho*nu) + "mu_" + str(init_particles) + "particles" + "/durham_sim"
if not os.path.isdir(outpath+folderstr):
    os.makedirs(outpath+folderstr)
outfile_particles = pset.ParticleFile(name=outpath + folderstr + "/particles", outputdt=timedelta(seconds=dt_NS))#, type='indexed')

for steps in range(num_steps):
    t += dt_NS

    # Load next firedrake field
    chk_in = checkpointing.HDF5File(path_to_fields + "u_func_" + "%05.0f" % (steps+1) + ".h5", file_mode='r')
    chk_in.read(up_sol, "/up")
    chk_in.close()


    if (t >= t_add_particles):  # and (steps%parcels_interval==0):
        # Extract velocity field at particle positions
        for particle in pset:
            try:
                particle.u, particle.v, particle.w = u_sol.at([particle.lon, particle.lat, particle.depth])
            except PointNotInDomainError:
                tqdm.write("\nRemoved out-of-domain particle %d at point (%f, %f, %f) during velocity extraction." % (
                particle.id, particle.lon, particle.lat, particle.depth))
                particle.delete()

        # Update particle positions
        pset.execute(AdvectionEE_firedrake_3D,  # the kernel (which defines how particles move)
                     runtime=timedelta(seconds=dt_NS * parcels_interval),  # the total length of the run
                     dt=timedelta(seconds=dt_NS),  # the timestep of the kernel
                     recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle},
                     output_file=outfile_particles)


print("Ended")
