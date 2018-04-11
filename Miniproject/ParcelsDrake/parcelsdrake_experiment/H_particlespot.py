# Imports
from datetime import timedelta

import numpy as np
from firedrake import *
from parcels import *
from tqdm import tqdm

from incflow.inc_navier_stokes import IncNavierStokes
from pdist import uniform_H_dist

def extract_field(u, grid_incr):
    u_field_1 = np.array(
        [[u.at(x, y) for x in np.arange(0, 2.8, grid_incr)] for y in np.arange(0, 13.6, grid_incr)])
    u_field_2 = np.vstack(([[(0, 0) for x in np.arange(2.8 + grid_incr, 10.8, grid_incr)] for y in
                            np.arange(8.8, 13.6, grid_incr)],
                           [[u.at([x, y]) for x in np.arange(2.8 + grid_incr, 10.8, grid_incr)] for y in
                            np.arange(4.8, 8.8 - grid_incr, grid_incr)],
                           [[(0, 0) for x in np.arange(2.8 + grid_incr, 10.8, grid_incr)] for y in
                            np.arange(0, 4.8 - grid_incr, grid_incr)]))
    u_field_3 = np.array([[u.at(x, y) for x in np.arange(10.8 + grid_incr, 13.6, grid_incr)] for y in
                          np.arange(0, 13.6, grid_incr)])
    u_field = np.hstack((u_field_1, u_field_2, u_field_3))
    return u_field

# Navier-Stokes firedrake setup
mesh_NS = Mesh("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/meshes/H.msh")
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
#up_init_path = "/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/inputs/60s_0.001mu_H.h5"
INS.setup_solver()#up_init=up_init_path)
u_sol, p_sol = INS.up.split()


# Lagrangian parcels setup
Bval = 0  # gyrotactic reorientation timescale
num_particles = 500
init_particles = 100
max_replace = 50
repeat_particles = 20
grid_res = 100
grid_start = 0
grid_end = 13.6
grid_incr = abs(grid_start - grid_end)/grid_res
lon = np.arange(grid_start, grid_end, grid_incr)
lat = np.arange(grid_start, grid_end, grid_incr)
grid = RectilinearZGrid(lon=lon, lat=lat, time=np.zeros(1), mesh='flat')

Ufield_init = Field(name='U', data=np.zeros((len(lon), len(lat))), time=0, grid=grid, allow_time_extrapolation=True)
Vfield_init = Field(name='V', data=np.zeros((len(lon), len(lat))), time=0, grid=grid, allow_time_extrapolation=True)
fieldset = FieldSet(U=Ufield_init, V=Vfield_init)


class FiredrakeParticle(JITParticle):
    u = Variable('u', dtype=np.float32, initial=0.)
    v = Variable('v', dtype=np.float32, initial=0.)


pfield_init = Field(name='pfield_init', data=uniform_H_dist(grid), grid=grid)
pfield_replace = Field(name='pfield_replace', data=uniform_H_dist(grid), grid=grid)
# pfield_replace = Field(name="replacefield", data=np.vstack((np.concatenate((np.zeros(20), np.repeat(0.1, 3), np.zeros(54), np.repeat(0.1, 2), np.zeros(21))),
#                                                             np.zeros((99, 100)))),
#                      grid=grid)
repeatdt = 10
pset = ParticleSet(fieldset=fieldset, pclass=FiredrakeParticle, lon=np.repeat(2.8, repeat_particles), lat=np.linspace(4.8, 8.8, repeat_particles))#, repeatdt=repeatdt)
# pset = ParticleSet.from_line(fieldset=fieldset,
#                              pclass=FiredrakeParticle,
#                              start=(2.8, 4.8),
#                              finish=(2.8, 8.8),
#                              size=init_particles)
# pset = ParticleSet.from_field(fieldset=fieldset,
#                               pclass=FiredrakeParticle,
#                               start_field=pfield_init,
#                               size=init_particles)

# Custom Parcels kernels
def DeleteParticle(particle, fieldset, time, dt):  # delete particles who run out of bounds.
    particle.delete()


# Time loop setup
t = 0
t_add_particles = 60
t_end = 120
parcels_interval = 1
num_steps = int((t_end - t)/INS.dt)
outfile_u = File("/media/alexander/DATA/Ubuntu/Miniproject/Firedrake Learning/outputs/NS_H/Julian/0.0Ti_60Tf_0.1dt_0.001mu_velBC/postparticles/u.pvd")
for steps in tqdm(range(num_steps)):
    t += dt_NS

    # Update firedrake solution
    u_sol, p_sol = INS.step()
    u_vort = curl(u_sol)

    # Update parcels field for plotting (slow so dont do often!).
    if t>=t_add_particles and steps%10 == 0:
        field_data = extract_field(u_sol, grid_incr=grid_incr)
        plotfieldset = FieldSet.from_data(data={'U': field_data[:, :, 0], 'V': field_data[:, :, 1]},
                                            dimensions={'lon': np.arange(grid_start, grid_end-grid_incr, grid_incr), 'lat': np.arange(grid_start, grid_end, grid_incr), 'time': np.array([t], float)},
                                            transpose=False, mesh='flat')

    if (t>=t_add_particles): #and (steps%parcels_interval==0):
        # Replace lost particles
        # if pset.size < num_particles:
        #     pset_replace = ParticleSet.from_line(fieldset=fieldset,
        #                                            pclass=FiredrakeParticle,
        #                                            start=(2.8, 4.8),
        #                                            finish=(2.8, 8.8),
        #                                            size=abs(min(num_particles - pset.size, max_replace)))
        #
        #     # pset_replace = ParticleSet.from_field(fieldset=fieldset,
        #     #                                       pclass=FiredrakeParticle,
        #     #                                       start_field=pfield_replace,
        #     #                                       size=abs(min(num_particles - pset.size, max_replace)))
        #     pset.add(pset_replace)

        # Extract velocity field at particle positions
        for particle in pset:
            try:
                particle.u, particle.v = u_sol.at([particle.lon, particle.lat])
                particle.vort_i, particle.vort_j = u_vort.at([particle.lon, particle.lat])
            except PointNotInDomainError:
                tqdm.write("\nRemoved out-of-domain particle %d at point (%f, %f) during velocity extraction." % (particle.id, particle.lon, particle.lat))
                particle.delete()
                continue

        # Update particle positions
        pset.execute(AdvectionEE_firedrake,  # the kernel (which defines how particles move)
                     runtime=timedelta(seconds=dt_NS * parcels_interval),  # the total length of the run
                     dt=timedelta(seconds=dt_NS),  # the timestep of the kernel
                     recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})
                     #moviedt=timedelta(seconds=dt_NS),
                     #movie_background_field=fieldset.V)

        # Save image to disk
        outfile_u.write(u_sol)
        pset.show(
            #savefile='/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/outputs/H/particlespot/animations/turbulencetest/img_debug/REPturbulencetest_' + str(steps).zfill(4),
            savefile="/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/outputs/H/particlespot/animations/sidebyside/img_particles/particles"+str(steps).zfill(4),
            field=plotfieldset.U,
            domain=(13.6, 0., 13.6, 0.),
            vmin=-1.1,
            vmax=1.1)

# Checkpoint final state.
chk_out = checkpointing.HDF5File("/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/outputs/H/fields/" + str(t_end) + "s_" + str(rho*nu) + "mu_" + "H" + ".h5", file_mode='w')
chk_out.write(u_sol, "/velocity")
chk_out.write(p_sol, "/pressure")
chk_out.write(INS.up, "/up")
chk_out.close()

print("Ended")