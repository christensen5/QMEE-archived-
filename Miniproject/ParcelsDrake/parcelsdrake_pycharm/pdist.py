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