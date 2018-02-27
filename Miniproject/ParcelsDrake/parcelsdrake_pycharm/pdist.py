import numpy as np


def uniform_H_Dist(grid):
    """Produce a uniform distribution over the entire H (zero values on rest of grid) with specified density.

    :param grid: Grid on which to define distribution.

    :returns dist_array: numpy array containing the distribution
    """

    if not grid.xdim==grid.ydim:
        raise ValueError('Grid is not square. Non-square grids not yet supported!')

    grid_res = grid.xdim

    dist_array = np.ones((grid_res, grid_res))
    for row in range(grid_res):
        for col in range(grid_res):
            if (row <= (2.8/13.6)*grid_res or row >= (8.8/13.6)*grid_res) and ((4.8/13.6)*grid_res <= col <= (10.8/13.6)*grid_res):
                dist_array[row, col] = 0

    return dist_array