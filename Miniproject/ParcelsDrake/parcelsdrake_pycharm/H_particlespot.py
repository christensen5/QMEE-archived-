# Imports
from datetime import timedelta

import numpy as np
from firedrake import *
from parcels import *
from tqdm import tqdm

from incflow.inc_navier_stokes import IncNavierStokes
from pdist import uniform_H_Dist

def extract_field(u, grid_res):
    u_field_1 = np.array(
        [[u.at(x, y) for x in np.arange(0, 2.8, 13.6 / 100)] for y in np.arange(0, 13.6 - 13.6 / 100, 13.6 / 100)])
    u_field_2 = np.vstack(([[(0, 0) for x in np.arange(2.8 + 13.6 / 100, 10.8, 13.6 / 100)] for y in
                            np.arange(8.8, 13.6 - 13.6 / 100, 13.6 / 100)],
                           [[u.at([x, y]) for x in np.arange(2.8 + 13.6 / 100, 10.8, 13.6 / 100)] for y in
                            np.arange(4.8, 8.8 - 13.6 / 100, 13.6 / 100)],
                           [[(0, 0) for x in np.arange(2.8 + 13.6 / 100, 10.8, 13.6 / 100)] for y in
                            np.arange(0, 4.8 - 13.6 / 100, 13.6 / 100)]))
    u_field_3 = np.array([[u.at(x, y) for x in np.arange(10.8 + 13.6 / 100, 13.6, 13.6 / 100)] for y in
                          np.arange(0, 13.6 - 13.6 / 100, 13.6 / 100)])
    u_field = np.hstack((u_field_1, u_field_2, u_field_3))
    return u_field

# Navier-Stokes firedrake setup
mesh_NS = Mesh("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/meshes/H_low.msh")
rho = 1
nu = 0.001
dt_NS = 0.1
INS = IncNavierStokes(mesh_NS, nu, rho, dt_NS)
W = INS.get_mixed_fs()
u_bcs = [DirichletBC(W.sub(0), Constant((0.0, 0.0)), 15),
         DirichletBC(W.sub(0), Constant((0.0, 1.0)), 16),  # inflow velocity of (0,1)
         DirichletBC(W.sub(0), Constant((0.0, 1.0)), 17)]  # outflow velocity of (0,1)
p_bcs = [DirichletBC(W.sub(1), Constant(0.0), 17)]  # outflow pressure of 0
INS.set_bcs(u_bcs, p_bcs)
INS.setup_solver()
u_sol, p_sol = INS.up.split()


# Lagrangian parcels setup
num_particles = 100
init_particles = 500
max_replace = 100
grid_res = 100
grid_start = 0
grid_end = 13.6
grid_incr = abs(grid_start - grid_end)/grid_res
lon = np.arange(grid_start, grid_end+grid_incr, grid_incr)
lat = np.arange(grid_start, grid_end+grid_incr, grid_incr)
grid = RectilinearZGrid(lon=lon, lat=lat, time=np.zeros(1), mesh='flat')

Ufield_init = Field(name='U', data=np.zeros((len(lon), len(lat))), time=0, grid=grid, allow_time_extrapolation=True)
Vfield_init = Field(name='V', data=np.zeros((len(lon), len(lat))), time=0, grid=grid, allow_time_extrapolation=True)
fieldset = FieldSet(U=Ufield_init, V=Vfield_init)


class FiredrakeParticle(ScipyParticle):
    u = Variable('u', dtype=np.float32, initial=0.)
    v = Variable('v', dtype=np.float32, initial=0.)


pfield_init = Field(name='pfield_init', data=uniform_H_Dist(grid), grid=grid)
pfield_replace = Field(name='pfield_replace', data=uniform_H_Dist(grid), grid=grid)
# pfield_replace = Field(name="replacefield", data=np.vstack((np.concatenate((np.zeros(20), np.repeat(0.1, 3), np.zeros(54), np.repeat(0.1, 2), np.zeros(21))),
#                                                             np.zeros((99, 100)))),
#                      grid=grid)#

pset = ParticleSet.from_field(fieldset=fieldset,
                              pclass=FiredrakeParticle,
                              start_field=pfield_init,
                              size=init_particles)

# Custom Parcels kernels
def DeleteParticle(particle, fieldset, time, dt):  # delete particles who run out of bounds.
    particle.delete()


# def AdvectionEE_firedrake(particle, fieldset, time, dt):
#     """Advection of particles using Explicit Euler (aka Euler Forward) integration,
#     and with particle velocities already provided from a Firedrake field."""
#     if (particle.lon < grid_start or particle.lat < grid_start or particle.lon > grid_end or particle.lat > grid_end) or (2.8 < particle.lon < 10.8 and (particle.lat < 4.8 or particle.lat > 8.8)):
#         particle.delete()
#         print("Deleted out-of-bounds particle %d at line 74." % particle.id)
#     else:
#         particle.lon += particle.u * dt
#         particle.lat += particle.v * dt


# Time loop setup
t = 0
t_add_particles = 60
t_end = 120  # 120s vortex form + 600s particles
parcels_interval = 1
num_steps = int((t_end - t)/INS.dt)

for steps in tqdm(range(num_steps)):
    t += dt_NS

    # noinspection PyRedeclaration
    u_sol, p_sol = INS.step()

    if (t>=t_add_particles) and (steps%parcels_interval==0):
        # Replace lost particles
        if pset.size < num_particles:
            pset_replace = ParticleSet.from_field(fieldset=fieldset,
                                                  pclass=FiredrakeParticle,
                                                  start_field=pfield_replace,
                                                  size=abs(min(num_particles - pset.size, max_replace)),
                                                  time=t)
            pset.add(pset_replace)

        # Extract particle grid positions & update Parcels' velocity field
        # Uarray = fieldset.U.data
        # Varray = fieldset.V.data
        # for particle in pset.particles:
        #     grid_x, grid_y = [(particle.lon/(13.6/grid_res)), (particle.lat/(13.6/grid_res))]
        #     try:
        #         particle_u, particle_v = u_sol.at([particle.lon, particle.lat])
        #     except PointNotInDomainError:
        #         print("\nDeleted out-of-domain particle %d from point (%f, %f)." %(particle.id, particle.lon, particle.lat))
        #         particle.delete()
        #         continue
        #     Uarray[0, floor(grid_y):ceil(grid_y), floor(grid_x):ceil(grid_x)] = particle_u
        #     Varray[0, floor(grid_y):ceil(grid_y), floor(grid_x):ceil(grid_x)] = particle_v
        #particle_x = np.zeros([])
        for particle in range(pset.size):
            try:
                pset[particle].u, pset[particle].v = u_sol.at([pset[particle].lon, pset[particle].lat])
            except PointNotInDomainError:
                tqdm.write("\nRemoved out-of-domain particle %d at point (%f, %f) during velocity extraction." % (pset[particle].id, pset[particle].lon, pset[particle].lat))
                pset[particle].delete()
                continue
            #particle_x = np.array([pset[particle].lon, pset[particle].lat]) if particle_x.shape==() else np.vstack((particle_x, [pset[particle].lon, pset[particle].lat]))

        #grid = RectilinearZGrid(lon=lon, lat=lat, time=np.array([t]), mesh='flat')
        #Ufield = Field(name='U', data=Uarray, grid=grid, allow_time_extrapolation=True)
        #Vfield = Field(name='V', data=Varray, grid=grid, allow_time_extrapolation=True)
        #fieldset.advancetime(FieldSet(U=Ufield, V=Vfield))
        pset.execute(AdvectionEE_firedrake,  # the kernel (which defines how particles move)
                     runtime=timedelta(seconds=dt_NS * parcels_interval),  # the total length of the run
                     dt=timedelta(seconds=dt_NS),  # the timestep of the kernel
                     recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})
                     #moviedt=timedelta(seconds=dt_NS),
                     #movie_background_field=fieldset.V)

        # Save image to disk
        # pset.show(
        #     savefile='/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/outputs/H/particlespot/animations/kerneltest_' +
        #              str(steps).zfill(4), field=fieldset.U, domain=(13.6, 0., 13.6, 0.), vmin=-1.1, vmax=1.1)
        # pset.show(savefile='/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/outputs/H/comparison/spot/spot' + str(steps).zfill(4), field=fieldset.U, vmin=-1.1, vmax=1.1)

print("Ended")