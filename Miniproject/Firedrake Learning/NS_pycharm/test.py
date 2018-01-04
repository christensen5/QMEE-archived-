from firedrake import *
import matplotlib.pyplot as plt


# mesh = UnitSquareMesh(32, 32)
# BDM = FunctionSpace(mesh, "BDM", 1)
# DG = FunctionSpace(mesh, "DG", 0)
# W = BDM * DG
# sigma, u = TrialFunctions(W)
# tau, v = TestFunctions(W)
# x, y = SpatialCoordinate(mesh)
# f = Function(DG).interpolate(
#     10*exp(-(pow(x - 0.5, 2) + pow(y - 0.5, 2)) / 0.02))
#
# a = (dot(sigma, tau) + div(tau)*u + div(sigma)*v)*dx
# L = - f*v*dx

mesh = UnitSquareMesh(16, 16)
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)
Z = V * Q

u = TrialFunction(V)
p = TrialFunction(Q)
v = TestFunction(V)
q = TestFunction(Q)

u_star = Function(V)
u_now = Function(V)
f = Constant(as_vector([0.0, 0.0]))

F = inner((u - u_now), v) * dx - inner(f, v) * dx
a, L = system(F)

solve(a == L, u_star) #, bcs=bcu)

