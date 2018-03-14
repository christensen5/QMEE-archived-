import netCDF4
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

nc = netCDF4.Dataset("/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/outputs/H_3D/particlespot/0Ti_2Tf_0.1dt_0.01mu/3D_test.nc")

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_xlabel("Longitude")
ax.set_ylabel("Depth")
ax.set_zlabel("Latitude")

for particle in range(nc.variables['time'].shape[0]):
    x = nc.variables["lon"][particle]
    y = -nc.variables["z"][particle]
    z = nc.variables["lat"][particle]

    #plot = ax.scatter(x, y, z, c=z, s=20, marker="o")
    plot = ax.plot(x, y, z, 'o-')

plt.show()
