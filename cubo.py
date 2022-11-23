# Preliminaries and mesh
from dolfin import *
import numpy.linalg as LA
#from mshr import *
from fenics import *
import sympy as sp
import numpy as np
from mpi4py import MPI as _MPI

mesh = Mesh('cubo.xml') 


mesh_f = File ("./new/cubo.pvd")
mesh_f << mesh

# Initialization of the iterative procedure and output requests
deltaT  = 1.0
kappa=1E-5

# Define Space
V = FunctionSpace(mesh, 'CG', 1)
W = VectorFunctionSpace(mesh, 'CG', 1)
WW = FunctionSpace(mesh, 'DG', 0)
p, q = TrialFunction(V), TestFunction(V)
pnew, pold, Hold = Function(V), Function(V), Function(V)

# ############ E ##############################

T = FunctionSpace(mesh, 'CG', 1)

local_range = T.dofmap().ownership_range()
local_dim = local_range[1] - local_range[0]
file_vals = open ( 'E_cubo.txt' )
list_vals = [ float ( line ) for line in file_vals.readlines() ]

array_vals = np.array ( list_vals)
E = Function ( T )
d2v = dof_to_vertex_map ( T )
global_vertex_numbers = mesh.topology().global_indices(0)
global_vertices = global_vertex_numbers[d2v[:local_dim]]
local_data = array_vals[global_vertices]
E.vector()[:] = local_data[:local_dim]
file = File ("./0/E.pvd")
file << E

############ NU ##############################

TT = FunctionSpace(mesh, 'CG', 1)

local_range_nuu = TT.dofmap().ownership_range()
local_dim_nuu = local_range_nuu[1] - local_range_nuu[0]
file_vals_nuu = open ( 'nuu_cubo.txt' )
list_vals_nuu = [ float ( line ) for line in file_vals_nuu.readlines() ]
array_vals_nuu = np.array ( list_vals_nuu)

nuu = Function ( TT )
d2v_nuu = dof_to_vertex_map ( TT )
global_vertex_numbers_nuu = mesh.topology().global_indices(0)
global_vertices_nuu = global_vertex_numbers_nuu[d2v_nuu[:local_dim_nuu]]
local_data_nuu = array_vals[global_vertices_nuu]
nuu.vector()[:] = local_data_nuu[:local_dim_nuu]
file = File ("./1/nuu.pvd")
file << nuu



