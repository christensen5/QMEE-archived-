from firedrake import *
import numpy as np

mesh = Mesh("/home/alexander/Documents/QMEE/Miniproject/Firedrake Learning/meshes/H_low.msh")
V = VectorFunctionSpace(mesh, "CG", 2)

u = Function(V)
chk_in = checkpointing.HDF5File(
    "/media/alexander/DATA/Ubuntu/Miniproject/Firedrake Learning/outputs/NS_H/Julian/60T_0.05dt_0.01mu/final.h5",
    file_mode='r')
chk_in.read(u, "/velocity")
chk_in.close()

# Generate cartesian velocity field
u_field_1 = np.array([[u.at(x, y) for x in np.arange(0, 2.8, 13.6/100)] for y in np.arange(0, 13.6 - 13.6/100, 13.6/100)])
u_field_2 = np.vstack(([[(0, 0) for x in np.arange(2.8 + 13.6/100, 10.8, 13.6/100)] for y in np.arange(8.8, 13.6 - 13.6/100, 13.6/100)],
                      [[u.at([x, y]) for x in np.arange(2.8 + 13.6/100, 10.8, 13.6/100)] for y in np.arange(4.8, 8.8 - 13.6/100, 13.6/100)],
                      [[(0, 0) for x in np.arange(2.8 + 13.6/100, 10.8, 13.6/100)] for y in np.arange(0, 4.8 - 13.6/100, 13.6/100)]))
u_field_3 = np.array([[u.at(x, y) for x in np.arange(10.8 + 13.6/100, 13.6, 13.6/100)] for y in np.arange(0, 13.6 - 13.6/100, 13.6/100)])
u_field = np.hstack((u_field_1, u_field_2, u_field_3))

np.save("/media/alexander/DATA/Ubuntu/Miniproject/Parcels Learning/fields/offline_grid_H_low.npy", u_field)