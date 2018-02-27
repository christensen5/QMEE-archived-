# Imports
from firedrake import *
from incflow import *
from incflow.inc_navier_stokes import IncNavierStokes
from parcels import *
import os, shutil
from tqdm import tqdm
import numpy as np
from math import floor, ceil
from datetime import timedelta
from pdist import uniform_H_Dist


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


# Lagrangian parcels setup
num_particles = 10
init_particles = 5
grid_res = 100
lon = np.arange(0, 13.6, 13.6/grid_res)
lat = np.arange(0, 13.6, 13.6/grid_res)
grid = RectilinearZGrid(lon=lon, lat=lat, time=np.zeros(1), mesh='flat')

Ufield_init = Field(name='U', data=np.zeros((len(lon), len(lat))), time=0, grid=grid, allow_time_extrapolation=True)
Vfield_init = Field(name='V', data=np.zeros((len(lon), len(lat))), time=0, grid=grid, allow_time_extrapolation=True)
fieldset = FieldSet(U=Ufield_init, V=Vfield_init)

pfield_init = Field(name='pfield_init',
                             data=uniform_H_Dist(grid),
                             grid=grid)
pfield_replace = Field(name='pfield_replace',
                                data=uniform_H_Dist(grid),
                                grid=grid)
pset = ParticleSet.from_field(fieldset=fieldset,
                              pclass=JITParticle,
                              start_field=pfield_init,
                              size=init_particles)

# Custom Parcels kernels
def DeleteParticle(particle, fieldset, time, dt):  # delete particles who run out of bounds.
    particle.delete()


# Time loop setup
step = 0
t = 0
t_add_particles = 1
t_end = 2  # 120s vortex form + 600s particles
parcels_interval = 1
num_steps = int((t_end - t)/INS.dt)

for steps in tqdm(range(num_steps)):
    t += dt_NS

    u_sol, p_sol = INS.step()

    if (t>=t_add_particles) and (steps%parcels_interval==0):
        # Replace lost particles
        if pset.size < num_particles:
            pset_replace = ParticleSet.from_field(fieldset=fieldset,
                                                  pclass=JITParticle,
                                                  start_field=pfield_replace,
                                                  size=min(num_particles - pset.size, 1))
            pset.add(pset_replace)

        # Save image to disk
        #pset.show(savefile='/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/outputs/H/particlespot/animations/'+str(steps).zfill(4), field=fieldset.U, vmin=-1.1, vmax=1.1)
        #pset.show(savefile='/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/outputs/H/comparison/spot/spot' + str(steps).zfill(4), field=fieldset.U, vmin=-1.1, vmax=1.1)


        # Extract particle grid positions & update Parcels' velocity field
        Uarray = fieldset.U.data
        Varray = fieldset.V.data
        for particle in pset.particles:
            grid_x, grid_y = [(particle.lon/(13.6/grid_res)), (particle.lat/(13.6/grid_res))]
            try:
                particle_u, particle_v = u_sol.at([particle.lon, particle.lat])
            except PointNotInDomainError:
                print("\nDeleted particle %d from point (%f, %f)." %(particle.id, particle.lon, particle.lat))
                particle.delete()
                continue
            Uarray[0, floor(grid_y):ceil(grid_y), floor(grid_x):ceil(grid_x)] = particle_u
            Varray[0, floor(grid_y):ceil(grid_y), floor(grid_x):ceil(grid_x)] = particle_v
        grid = RectilinearZGrid(lon=lon, lat=lat, time=np.array([t]), mesh='flat')
        Ufield = Field(name='U', data=Uarray, grid=grid, allow_time_extrapolation=True)
        Vfield = Field(name='V', data=Varray, grid=grid, allow_time_extrapolation=True)
        fieldset.advancetime(FieldSet(U=Ufield, V=Vfield))

        pset.execute(AdvectionRK4,  # the kernel (which defines how particles move)
                     runtime=timedelta(seconds=dt_NS * parcels_interval),  # the total length of the run
                     dt=timedelta(seconds=dt_NS),  # the timestep of the kernel
                     recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})
                     #moviedt=timedelta(seconds=dt_NS),
                     #movie_background_field=fieldset.V)

print("Ended")