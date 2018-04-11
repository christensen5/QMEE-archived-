# Imports
from datetime import timedelta
import os

import numpy as np
from firedrake import *
from parcels import *
from tqdm import tqdm

from incflow.inc_navier_stokes_3D import IncNavierStokes3D
from pdist import uniform_H_dist, uniform_H_bar_dist, uniform_H_bar_entry_dist

def extract_field_3D(u, grid_start, grid_end, grid_incr, depth):
    """ Produce a discrete velocity field (discretised according to the provided grid) from firedrake field u """
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
#mesh_NS = Mesh("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/meshes/H_3D_high.msh")
rho = 1
nu = 0.01
dt_NS = 0.5
INS = IncNavierStokes3D(mesh_NS, nu, rho, dt_NS)
W = INS.get_mixed_fs()
up_sol = Function(W)
path_to_fields = "/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/inputs/Hmixed/data1/0.5dt_0.01mu/h5/"
#path_to_fields = "/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/inputs/Hhigh/data1h/0.5dt_0.01mu/h5/"
chk_in = checkpointing.HDF5File(path_to_fields + "u_func_00001.h5", file_mode='r')
chk_in.read(up_sol, "/up")
chk_in.close()
u_sol, p_sol = up_sol.split()


# Lagrangian parcels setup
B_val = float(2.0)
Swim_val = 1  # REMOVE ONCE RANDOMISED WITHIN FIREDRAKEPARTICLE CLASS
num_particles = 500
init_particles = 2000
max_replace = 50
repeat_particles = 100
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


class MotileParticle3D(JITParticle):
    u = Variable('u', dtype=np.float32, initial=0.)
    v = Variable('v', dtype=np.float32, initial=0.)
    w = Variable('w', dtype=np.float32, initial=0.)
    vort_i = Variable('vort_i', dtype=np.float32, initial=0.)
    vort_j = Variable('vort_j', dtype=np.float32, initial=0.)
    vort_k = Variable('vort_k', dtype=np.float32, initial=0.)
    dir_x = Variable('dir_x', dtype=np.float32, initial=0.)  # RANDOMISE INITIAL VALUE
    dir_y = Variable('dir_y', dtype=np.float32, initial=0.)  # RANDOMISE INITIAL VALUE
    dir_z = Variable('dir_z', dtype=np.float32, initial=0.)  # RANDOMISE INITIAL VALUE
    B = Variable('B', dtype=np.float32, initial=B_val)
    v_swim = Variable('v_swim', dtype=np.float32, initial=Swim_val)  # RANDOMISE INITIAL VALUE


# pfield_uniform = Field(name='pfield_uniform', data=uniform_H_dist(grid), transpose=False, grid=grid)
# pset = ParticleSet.from_field(fieldset=fieldset,
#                               pclass=FiredrakeParticle3D,
#                               start_field=pfield_uniform,
#                               depth=np.random.rand(init_particles) * depth,
#                               size=init_particles)
pfield_Hbar_uniform = Field(name='pfield_Hbar_uniform', data=uniform_H_bar_dist(grid), transpose=False, grid=grid)
pset = ParticleSet.from_field(fieldset=fieldset,
                              pclass=FiredrakeParticle3D,
                              start_field=pfield_Hbar_uniform,
                              depth=np.random.rand(init_particles) * depth,
                              size=init_particles)
# pfield_Hbar_entry = Field(name='pfield_Hbar_entry', data=uniform_H_bar_entry_dist(grid), transpose=False, grid=grid)
# pset = ParticleSet.from_field(fieldset=fieldset,
#                               pclass=FiredrakeParticle3D,
#                               start_field=pfield_Hbar_entry,
#                               depth=np.linspace(0, depth, init_particles),
#                               size=init_particles)


# Custom Parcels kernels
def DeleteParticle(particle, fieldset, time, dt):  # delete particles who run out of bounds.
    particle.delete()


def Gyrotaxis(particle, fieldset, time, dt):  # Gyrotactically-determined alignment and swimming of particles.
    (u1, v1) = fieldset.UV[time, particle.lon, particle.lat, particle.depth] # ad-hoc check if particle has left the square - works fine BUT sub-optimal

    # # Re-align the particle
    di = 0.5 * ((1/particle.B * -1 * particle.dir_z * particle.dir_x) + (particle.vort_j * particle.dir_z) - (particle.vort_k * particle.dir_y))
    dj = 0.5 * ((1/particle.B * -1 * particle.dir_z * particle.dir_y) + (particle.vort_k * particle.dir_x) - (particle.vort_i * particle.dir_z))
    dk = 0.5 * ((1/particle.B * (1 - particle.dir_z**2)) + (particle.vort_i * particle.dir_y) - (particle.vort_j * particle.dir_x))
    newdir = np.array([particle.dir_x + di, particle.dir_y + dj, particle.dir_z + dk])

    particle.dir_x, particle.dir_y, particle.dir_z = newdir / (newdir[0]**2 + newdir[1]**2 + newdir[2]**2)**0.5  # alignment vector must be unit-length

    # Update position (FIELD + SWIM)
    particle.lon += (particle.u * dt) + (particle.dir_x * particle.v_swim * dt)
    particle.lat += (particle.v * dt) + (particle.dir_y * particle.v_swim * dt)
    particle.depth += (particle.w * dt) + (particle.dir_z * particle.v_swim * dt)


# Time loop setup
removals = 0
removals_low = 0
removed = set()
t = 0
t_add_particles = 300
t_end = 420
parcels_interval = 1
num_steps = int((t_end - t)/dt_NS)
outpath = "/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/outputs/H_3D/Hmixed/particlespot/offline/"
#outpath = "/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/outputs/H_3D/Hhigh/particlespot/offline/"
folderstr = str(t) + "Ti_" + str(t_end) + "Tf_" + str(dt_NS) + "dt_" + str(rho*nu) + "mu_" + str(init_particles) + "particles" + "/durham_sim"
if not os.path.isdir(outpath+folderstr):
    os.makedirs(outpath+folderstr)
outfile_particles = pset.ParticleFile(name=outpath + folderstr + "/particles", outputdt=timedelta(seconds=dt_NS))#, type='indexed')
for steps in tqdm(range(num_steps)):
    t += dt_NS

    # Load next firedrake field
    chk_in = checkpointing.HDF5File(path_to_fields + "u_func_" + "%05.0f" % (steps + 1) + ".h5", file_mode='r')
    chk_in.read(up_sol, "/up")
    chk_in.close()


    if (t >= t_add_particles):  # and (steps%parcels_interval==0):
        # Extract velocity field at particle positions
        for particle in pset:
            try:
                particle.u, particle.v, particle.w = u_sol.at([particle.lon, particle.lat, particle.depth])
                #particle.vort_i, particle.vort_j, particle.vort_k = u_vort.at([particle.lon, particle.lat, particle.depth])
            except PointNotInDomainError:
                if particle.depth < 0:
                    particle.depth = 0.001
                if particle.depth > depth:
                    particle.depth = depth - 0.001
                if 2.8 < particle.lon < 10.8 and particle.lat > 8.8:
                    particle.lat = 8.8 - 0.001
                    particle.u, particle.v, particle.w = u_sol.at([particle.lon, particle.lat, particle.depth])
                    # particle.vort_i, particle.vort_j, particle.vort_k = u_vort.at([particle.lon, particle.lat, particle.depth])

                else:
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