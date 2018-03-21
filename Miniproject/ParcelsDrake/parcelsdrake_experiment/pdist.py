import numpy as np


def uniform_H_Dist(grid):
    """Produce a uniform distribution over the entire H (zero values on rest of grid) with specified density.

    :param grid: Grid on which to define distribution.

    :returns dist_array: numpy array containing the distribution
    """

    if not grid.xdim==grid.ydim:
        raise ValueError('Grid is not square. Non-square grids not yet supported!')

    grid_dim = grid.xdim

    dist_array = np.zeros((grid_dim, grid_dim))
    for row in range(grid_dim):
        for col in range(grid_dim):
            if col < (2.8/13.6)*grid_dim or col > (10.8/13.6)*grid_dim:
                dist_array[row, col] = 1
            elif ((4.8/13.6)*grid_dim < row < (8.8/13.6)*grid_dim):
                dist_array[row, col] = 1

    return dist_array


## CURRENTLY DEFUNCT - parcels initialised particleset_from_field in 3D by taking a 2D probability field and SEPARATE depth values.
# def uniform_H_Dist_3D(grid):
#     """Produce a uniform distribution over the entire 3D H (zero values on rest of grid) with specified density.
#
#         :param grid: Grid on which to define distribution.
#
#         :returns dist_array: numpy array containing the distribution
#         """
#
#     if not grid.xdim==grid.ydim:
#         raise ValueError('Grid is not square (w.r.t lon and lat). Non-square grids not yet supported!')
#
#     xdim = grid.xdim
#     ydim = grid.ydim
#     zdim = grid.zdim
#
#     dist_array = np.zeros((xdim, ydim, zdim))
#
#     for row in range(xdim):
#         for col in range(ydim):
#             if col < (2.8/13.6)*ydim or col > (10.8/13.6)*ydim:
#                 dist_array[row, col, :] = np.ones(zdim)
#             elif ((4.8/13.6)*xdim < row < (8.8/13.6)*xdim):
#                 dist_array[row, col, :] = np.ones(zdim)
#
#     return dist_array
