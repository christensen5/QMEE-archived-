# Imports
from firedrake import *
from incflow import *
from incflow.inc_navier_stokes import IncNavierStokes
from parcels import *
import os, shutil
from tqdm import tqdm
import numpy as np
from datetime import timedelta

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

def extract_field(u):
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


# Lagrangian parcels setup
num_particles = 10
init_particles = 5
lon = np.arange(0, 13.6-13.6/100, 13.6/100)
lat = np.arange(0, 13.6-13.6/100, 13.6/100)
u, p = INS.up.split()
field_data = extract_field(u)
fieldset = FieldSet.from_data(data={'U': field_data[:, :, 0], 'V': field_data[:, :, 1]}, dimensions={'lon': lon, 'lat': lat, 'time': np.array([0.0], float)}, transpose=False, mesh='flat')
# pset = ParticleSet.from_list(fieldset=fieldset,
#                              pclass=JITParticle,
#                              lon=[2.6, 2.6, 2.6, 2.7, 2.7, 2.7, 10.9, 10.9, 10.9, 11.0, 11.0, 11.0],
#                              lat=[1.5, 2.5, 3.5, 1.5, 2.5, 3.5, 1.5,  2.5,  3.5,  1.5,  2.5,  3.5],
#                              time=[10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0])
#                              # lon=np.hstack((np.arange(2.6, 2.8, 0.01), np.arange(10.81, 11.01, 0.01))).tolist(),
#                              # lat=np.full((40, 1), 0.5).tolist())
#pfield = Field(name='particlefield', data=np.load("/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/inputs/particle_H_central.npy"), lon=lon, lat=lat)
#pfield = Field(name="particlefield", data=np.vstack((np.concatenate(
#    (np.zeros(2), np.repeat(1.0 / 40, 19), np.zeros(57), np.repeat(1.0 / 40, 19), np.zeros(2))), np.zeros((98, 99)))),
#                     lon=lon, lat=lat)  # Uniform dist along bottom boundary
#replacefield = Field(name="replacefield", data=np.vstack((np.concatenate(
#    (np.zeros(20), np.repeat(0.1, 3), np.zeros(54), np.repeat(0.1, 2), np.zeros(20))), np.zeros((98, 99)))),
#                     lon=lon, lat=lat)
pfield = Field(name="particlefield", data=np.hstack((np.zeros((99, 22)),
                                                     np.vstack((np.zeros((36, 1)), np.repeat(1/34, 27).reshape(27, 1), np.zeros((36, 1)))),
                                                     np.zeros((99, 53)),
                                                     np.vstack((np.zeros((36, 1)), np.repeat(1/34, 27).reshape(27, 1), np.zeros((36, 1)))),
                                                     np.zeros((99, 22)))),
               lon=lon,
               lat=lat)
replacefield = Field(name="replacefield_bar_entry", data=np.hstack((np.zeros((99, 22)),
                                                     np.vstack((np.zeros((36, 1)), np.repeat(1/34, 27).reshape(27, 1), np.zeros((36, 1)))),
                                                     np.zeros((99, 53)),
                                                     np.vstack((np.zeros((36, 1)), np.repeat(1/34, 27).reshape(27, 1), np.zeros((36, 1)))),
                                                     np.zeros((99, 22)))),
                     lon=lon,
                     lat=lat)

pset = ParticleSet.from_field(fieldset=fieldset,
                              pclass=JITParticle,
                              start_field=pfield,
                              size=init_particles)


def DeleteParticle(particle, fieldset, time, dt):
    particle.delete()


# Time loop setup
step = 0
t = 0
t_add_particles = 60
t_end = 120  # 120s vortex form + 600s particles
parcels_interval = 1
num_steps = int((t_end - t)/INS.dt)

# simstr_NS = str(t_end) + "T_" + str(dt) + "dt_" + str(rho*nu) + "mu"
# folderstr = "/home/alexander/Documents/QMEE/Miniproject/ParcelsDrake/outputs/H/NS" + simstr_NS
# outfile_u = File(folderstr + "/u.pvd")
# outfile_p = File(folderstr + "/p.pvd")
#outfile_parcels = pset.ParticleFile(name="/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/outputs/H/testrun", outputdt=timedelta(seconds=dt_NS*10))

for steps in tqdm(range(num_steps)):
    t += INS.dt

    u_sol, p_sol = INS.step()

    if (t >= t_add_particles) and (steps%parcels_interval == 0):
        if pset.size < num_particles:
            pset_replace = ParticleSet.from_field(fieldset=fieldset,
                                                  pclass=JITParticle,
                                                  start_field=replacefield,
                                                  size=min(num_particles-pset.size, 1))
            pset.add(pset_replace)

        #pset.show(savefile='/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/outputs/H/wholefield/animations/vortex_entry/img_barEntry_10000ps/particles'+str(steps).zfill(4), field=fieldset.U, vmin=-1.1, vmax=1.1)
        pset.show(savefile='/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/outputs/H/comparison/whole/whole' + str(steps).zfill(4), field=fieldset.U, vmin=-1.1, vmax=1.1)

        field_data = extract_field(u_sol)
        fieldset.advancetime(FieldSet.from_data(data={'U': field_data[:, :, 0], 'V': field_data[:, :, 1]}, dimensions={'lon': lon, 'lat': lat, 'time': np.array([t], float)},
                           transpose=False, mesh='flat'))
        pset.execute(AdvectionRK4,  # the kernel (which defines how particles move)
                     runtime=timedelta(seconds=dt_NS * parcels_interval),  # the total length of the run
                     dt=timedelta(seconds=dt_NS),  # the timestep of the kernel
                     recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})
#                     moviedt=timedelta(seconds=dt_NS*parcels_interval),
#                     movie_background_field=fieldset.V)
#                     output_file=outfile_parcels)


# #    if steps % 10 == 0:
#     outfile_u.write(u_sol)
#     outfile_p.write(p_sol)

#plotTrajectoriesFile('../outputs/H/15T_500n_random_start.nc')
