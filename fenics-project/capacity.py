# -*- coding: utf-8 -*-
"""
FEniCS tutorial demo program: Poisson equation with Dirichlet conditions.
Test problem is chosen to give an exact solution at all nodes of the mesh.
  -Laplace(u) = f    in the unit square
            u = u_D  on the boundary
  u_D = 1 + x^2 + 2y^2
    f = -6
"""

from __future__ import print_function
from fenics import *
from mshr import *
import numpy as np

# Create mesh and define function space

square =  Rectangle(Point(0, 0), Point(1,1.5))
cylinder =  Circle(Point(0.5, 1.2), 0.1)
little_square = Rectangle(Point(0.35, 0.35), Point(0.65,0.65))
#domain = square - cylinder
domain = square - little_square - cylinder
mesh =  generate_mesh(domain, 128)
V = FunctionSpace(mesh, 'P', 1)
W = VectorFunctionSpace(mesh, 'P', 1)

test1 = False
capacity_01 = False
capacity_10 = True

if test1:
    u_D = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]', degree=2)


# Define boundary condition


tol = 1E-14

boundary_markers = FacetFunction('size_t', mesh)

class outer_0(SubDomain):
    def inside(self, x, on_boundary):
        return  near(x[0], 0, tol) or near(x[0],1, tol) \
            or near(x[1],0, tol) or near(x[1],1.5, tol)
outer0 = outer_0()
outer0.mark(boundary_markers, 0)
                                 
class little_square_1(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0]>0.3 and x[0]<0.7 and x[1]>0.3 and x[1]<0.7
l_square1 = little_square_1()
l_square1.mark(boundary_markers,1)

class cyl_2(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0]>0.1 and x[0]<0.7 and x[1]>1. and x[1]<1.4
cyl2 = cyl_2()
cyl2.mark(boundary_markers,2)

# the boundary surface elements 
ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)
             
if capacity_10:
    u_outer = Expression('0', degree=2)
    u_lsquare = Expression('1', degree=2)
    u_cyl = Expression('0', degree=2)

if capacity_01:                                 
    u_outer = Expression('0', degree=2)
    u_lsquare = Expression('0', degree=2)
    u_cyl = Expression('1', degree=2)

if capacity_01 or capacity_10:    
    bc_out = DirichletBC(V, u_outer, outer0)
    bc_ls = DirichletBC(V, u_lsquare, l_square1)
    bc_cyl = DirichletBC(V, u_cyl, cyl2)

if test1:
    bc_out = DirichletBC(V, u_D, outer0)
    bc_ls = DirichletBC(V, u_D, l_square1)
    bc_cyl = DirichletBC(V, u_D, cyl2) 

bcs = [bc_out, bc_ls, bc_cyl]

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)

if capacity_01 or capacity_10: 
    f = Constant(-0.0)

if test1:
    f = Constant(-6.0)
#f = Expression('exp(x[0]*x[0] + x[1]*x[1])', degree=2)
    
a = dot(grad(u), grad(v))*dx
L = f*v*dx

# Compute solution
u = Function(V)
solve(a == L, u, bcs)

# Plot solution and mesh
plot(u)
plot(mesh)

# Save solution to file in VTK format

if capacity_01 :
    vtkfile_phi = File('poisson/potential_01.pvd')
    vtkfile_e = File('poisson/field_01.pvd')

if capacity_10 :
    vtkfile_phi = File('poisson/potential_10.pvd')
    vtkfile_e = File('poisson/field_10.pvd')

if test1 :
    vtkfile_phi = File('poisson/potential_test.pvd')
    vtkfile_e = File('poisson/field_test.pvd')
    
vtkfile_phi << u

# Compute de electric field
grad_u = project(-grad(u),W)


vtkfile_e << grad_u

# Compute the charge integral

n = FacetNormal(mesh)
square_flux = -dot(grad(u),n)*ds(1)
Q_square = assemble(square_flux)

cyl_flux = -dot(grad(u),n)*ds(2)
Q_cyl = assemble(cyl_flux)

print('Charge_square = {0:.4f}, Charge_cyl = {1:.4f} ' .format(Q_square, Q_cyl))

if test1:
# Compute error in L2 norm
    error_L2 = errornorm(u_D, u, 'L2')

# Compute maximum error at vertices
    vertex_values_u_D = u_D.compute_vertex_values(mesh)
    vertex_values_u = u.compute_vertex_values(mesh)
    error_max = np.max(np.abs(vertex_values_u_D - vertex_values_u))

# Print errors
    print('error_L2  =', error_L2)
    print('error_max =', error_max)

# Hold plot
#interactive()
