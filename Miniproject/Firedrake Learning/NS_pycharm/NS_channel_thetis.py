from thetis import *

mesh = RectangleMesh(30, 10, 30, 10)  # elts then length

# Set up function space and function object to hold solution
P1_2d = FunctionSpace(mesh, 'CG', 2)
bathymetry_2d = Function(P1_2d, name='bathymetry')
bathymetry_2d.assign(10)

# Set solver options and create solver object
t_end = 50 # total duration in seconds
t_export = 0.1  # export interval in seconds

solver_obj = solver2d.FlowSolver2d(mesh, bathymetry_2d)
options = solver_obj.options
options.simulation_export_time = t_export
options.simulation_end_time = t_end
options.timestepper_type = 'CrankNicolson'
options.timestep = 0.1

# Name boundaries
left_bdary = 1
right_bdary = 2

# Define boundary conditions
bdary_flow = as_vector([2, 0])
bcs = {}
bcs[left_bdary] = {'uv': bdary_flow}
bcs[right_bdary] = {'uv': bdary_flow}

# Apply boundary conditions to solver object
solver_obj.bnd_functions['shallow_water'] = bcs

solver_obj.iterate()





