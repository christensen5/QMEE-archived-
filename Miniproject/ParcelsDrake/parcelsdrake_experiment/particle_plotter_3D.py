import netCDF4
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.patches import Rectangle, PathPatch
import matplotlib.pyplot as plt
import numpy.random as nprand
import numpy as np
from tqdm import tqdm

nc = netCDF4.Dataset("/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/outputs/H_3D/Hhigh/particlespot/offline/0Ti_960Tf_0.5dt_0.01mu_70000particles/durham_sim/particles.nc")

numparticles = 5000

def plot_trajectories():
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.set_xlabel("X - Longitude")
    ax.set_ylabel("Y - Depth")
    ax.set_zlabel("Z - Latitude")
    ax.set_xlim3d(0, 13.6)
    ax.set_ylim3d(0, 2.8)
    ax.set_zlim3d(0, 13.6)

    # Plot 3D H outline
    alpha = 0.5
    p1 = Rectangle([0, 0], 2.8, 13.6, alpha=alpha)
    p2 = Rectangle([2.8, 4.8], 8, 4, alpha=alpha)
    p3 = Rectangle([10.8, 0], 2.6, 13.6, alpha=alpha)
    p1_z = Rectangle([0, 0], 2.8, 13.6, alpha=alpha)
    p2_z = Rectangle([2.8, 4.8], 8, 4, alpha=alpha)
    p3_z = Rectangle([10.8, 0], 2.6, 13.6, alpha=alpha)
    p_top = Rectangle([2.8, 0], 8, 2.8, alpha=alpha)
    p_bottom = Rectangle([2.8, 0], 8, 2.8, alpha=alpha)
    ax.add_patch(p1)
    ax.add_patch(p2)
    ax.add_patch(p3)
    ax.add_patch(p1_z)
    ax.add_patch(p2_z)
    ax.add_patch(p3_z)
    ax.add_patch(p_top)
    ax.add_patch(p_bottom)
    art3d.pathpatch_2d_to_3d(p1, z=0, zdir="y")
    art3d.pathpatch_2d_to_3d(p2, z=0, zdir="y")
    art3d.pathpatch_2d_to_3d(p3, z=0, zdir="y")
    art3d.pathpatch_2d_to_3d(p1_z, z=2.8, zdir="y")
    art3d.pathpatch_2d_to_3d(p2_z, z=2.8, zdir="y")
    art3d.pathpatch_2d_to_3d(p3_z, z=2.8, zdir="y")
    art3d.pathpatch_2d_to_3d(p_top, z=4.8, zdir="z")
    art3d.pathpatch_2d_to_3d(p_bottom, z=8.8, zdir="z")

    random_particles = np.random.random_integers(0, numparticles - 1, 1999)

    for particle in random_particles:
        #if not np.ma.is_masked(nc.variables["lon"][particle]):
        x = nc.variables["lon"][particle]
        y = nc.variables["z"][particle]
        z = nc.variables["lat"][particle]

        plot = ax.plot(x, y, z, 'o-', linewidth=0.5, markersize=1.0)

    plt.show()


def plot_positions(particles, timestep):
    fig = plt.figure(figsize=plt.figaspect(0.5))

    # Front view
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    ax.set_title("Front View")
    ax.set_xlabel("X - Longitude")
    ax.set_ylabel("Y - Depth")
    ax.set_zlabel("Z - Latitude")
    ax.set_yticks([])
    ax.set_xlim(0, 13.6)
    ax.set_ylim(0, 2.8)
    ax.set_zlim(0, 13.6)
    ax.view_init(elev=5, azim=275)

    # Plot 3D H outline
    alpha = 0.2
    p1 = Rectangle([0, 0], 2.8, 13.6, alpha=alpha)
    p2 = Rectangle([2.8, 4.8], 8, 4, alpha=alpha)
    p3 = Rectangle([10.8, 0], 2.6, 13.6, alpha=alpha)
    p1_z = Rectangle([0, 0], 2.8, 13.6, alpha=alpha)
    p2_z = Rectangle([2.8, 4.8], 8, 4, alpha=alpha)
    p3_z = Rectangle([10.8, 0], 2.6, 13.6, alpha=alpha)
    p_top = Rectangle([2.8, 0], 8, 2.8, alpha=alpha)
    p_bottom = Rectangle([2.8, 0], 8, 2.8, alpha=alpha)
    ax.add_patch(p1)
    ax.add_patch(p2)
    ax.add_patch(p3)
    ax.add_patch(p1_z)
    ax.add_patch(p2_z)
    ax.add_patch(p3_z)
    ax.add_patch(p_top)
    ax.add_patch(p_bottom)
    art3d.pathpatch_2d_to_3d(p1, z=0, zdir="y")
    art3d.pathpatch_2d_to_3d(p2, z=0, zdir="y")
    art3d.pathpatch_2d_to_3d(p3, z=0, zdir="y")
    art3d.pathpatch_2d_to_3d(p1_z, z=2.8, zdir="y")
    art3d.pathpatch_2d_to_3d(p2_z, z=2.8, zdir="y")
    art3d.pathpatch_2d_to_3d(p3_z, z=2.8, zdir="y")
    art3d.pathpatch_2d_to_3d(p_top, z=4.8, zdir="z")
    art3d.pathpatch_2d_to_3d(p_bottom, z=8.8, zdir="z")

    # Plot particles
    for particle in particles:
        if not np.ma.is_masked(nc.variables["lon"][particle]):
            x = nc.variables["lon"][particle][timestep]
            y = nc.variables["z"][particle][timestep]
            z = nc.variables["lat"][particle][timestep]

            plot = ax.scatter(x, y, z, 'o-', c='r', s=0.3, alpha=0.1)

    # Side view
    ax = fig.add_subplot(1, 2, 2, projection='3d')
    ax.set_title("Side View")
    ax.set_xlabel("X - Longitude")
    ax.set_ylabel("Y - Depth")
    ax.set_zlabel("Z - Latitude")
    ax.set_xticks([])
    ax.set_xlim(0, 13.6)
    ax.set_ylim(0, 2.8)
    ax.set_zlim(0, 13.6)
    ax.view_init(elev=5, azim=5)

    # Plot 3D H outline
    alpha = 0.2
    p1 = Rectangle([0, 0], 2.8, 13.6, alpha=alpha)
    p2 = Rectangle([2.8, 4.8], 8, 4, alpha=alpha)
    p3 = Rectangle([10.8, 0], 2.6, 13.6, alpha=alpha)
    p1_z = Rectangle([0, 0], 2.8, 13.6, alpha=alpha)
    p2_z = Rectangle([2.8, 4.8], 8, 4, alpha=alpha)
    p3_z = Rectangle([10.8, 0], 2.6, 13.6, alpha=alpha)
    p_top = Rectangle([2.8, 0], 8, 2.8, alpha=alpha)
    p_bottom = Rectangle([2.8, 0], 8, 2.8, alpha=alpha)
    ax.add_patch(p1)
    ax.add_patch(p2)
    ax.add_patch(p3)
    ax.add_patch(p1_z)
    ax.add_patch(p2_z)
    ax.add_patch(p3_z)
    ax.add_patch(p_top)
    ax.add_patch(p_bottom)
    art3d.pathpatch_2d_to_3d(p1, z=0, zdir="y")
    art3d.pathpatch_2d_to_3d(p2, z=0, zdir="y")
    art3d.pathpatch_2d_to_3d(p3, z=0, zdir="y")
    art3d.pathpatch_2d_to_3d(p1_z, z=2.8, zdir="y")
    art3d.pathpatch_2d_to_3d(p2_z, z=2.8, zdir="y")
    art3d.pathpatch_2d_to_3d(p3_z, z=2.8, zdir="y")
    art3d.pathpatch_2d_to_3d(p_top, z=4.8, zdir="z")
    art3d.pathpatch_2d_to_3d(p_bottom, z=8.8, zdir="z")

    # Plot particles
    for particle in particles:
        if not np.ma.is_masked(nc.variables["lon"][particle]):
            x = nc.variables["lon"][particle][timestep]
            y = nc.variables["z"][particle][timestep]
            z = nc.variables["lat"][particle][timestep]

            plot = ax.scatter(x, y, z, 'o-', c='r', s=0.3, alpha=0.1)

    plt.savefig("/media/alexander/DATA/Ubuntu/Miniproject/ParcelsDrake/outputs/H_3D/Hhigh/particlespot/offline/0Ti_960Tf_0.5dt_0.01mu_70000particles/durham_sim/img/" + "%05.0f" % timestep + ".png")
    plt.close()


if __name__ == "__main__":
    #plot_trajectories()
    particles = np.random.random_integers(0, 70000 - 1, 5000)
    for timestep in tqdm(range(777)):
        plot_positions(particles, timestep)