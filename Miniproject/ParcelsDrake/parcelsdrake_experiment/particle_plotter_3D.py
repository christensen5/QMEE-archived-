import netCDF4
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy.random as nprand
import numpy as np

nc = netCDF4.Dataset("/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/outputs/H_3D/particlespot/offline/0Ti_600Tf_0.1dt_0.01mu/turbulence_test_Hbar_entry.nc")#first_particle_test.nc")

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_xlabel("Longitude")
ax.set_ylabel("Depth")
ax.set_zlabel("Latitude")

random_particles = range(500, 600)#np.random.random_integers(0, 1000, 1000)

for particle in random_particles:
    if not np.ma.is_masked(nc.variables["lon"][particle]):
        x = nc.variables["lon"][particle]
        y = -nc.variables["z"][particle]
        z = nc.variables["lat"][particle]

        #plot = ax.scatter(x, y, z, c=z, s=20, marker="o")
        plot = ax.plot(x, y, z, 'o-', linewidth=0.5, markersize=1.0)

plt.show()
