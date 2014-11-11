#!/usr/bin/env python

import os, sys

bemppInstallPrefix = "/home/mahir/"
sys.path.append(bemppInstallPrefix + "bempp/python")

from bempp.lib import *
import numpy as np
from cmath import exp
from math import sqrt

from bempp import single_layer_navier as sln
from bempp import double_layer_navier as dln
from bempp import simple_vector_spaces as svs

# Input data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Material parameters
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rho = 2000.0
mu = 10.0e6
nu = 0.25
lamb = 2*mu*nu/(1-2*nu)

# "Hysteretic material damping" - not implemented yet...
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
D = 0.02      # Material damping ratio
eta = 2*D     # "Hysteretic loss coefficient"

# frequency range
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
w = 10.0 * 2.0 * np.pi

# wave numbers
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
k_s = np.sqrt(rho*w**2/mu)
k_p = np.sqrt( (1-2*nu) / (2*(1-nu)) )*k_s

# Quadrature
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
accuracyOptions = createAccuracyOptions()
accuracyOptions.doubleRegular.setRelativeQuadratureOrder(1)
accuracyOptions.doubleSingular.setRelativeQuadratureOrder(1)
accuracyOptions.singleRegular.setRelativeQuadratureOrder(1)
quadStrategy = createNumericalQuadratureStrategy("float64", "complex128", accuracyOptions)
assemblyOptions = createAssemblyOptions()
context = createContext(quadStrategy, assemblyOptions)

# The boundary mesh
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
grid = createGridFactory().importGmshGrid("triangular", "square_1m.msh")
# grid = createGridFactory().importGmshGrid("triangular", "square_r10_b1_250mm.msh")

# The Dirichlet part of the boundary
segmentD = GridSegment.closedDomain(grid, 16)
# The Neumann part of the boundary
segmentN = segmentD.complement()

# Function spaces
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
spaceConstant = svs.createPiecewiseConstantVectorSpace(grid)
spaceDConstant = svs.createPiecewiseConstantVectorSpace(grid, segmentD)
spaceNConstant = svs.createPiecewiseConstantVectorSpace(grid, segmentN)
spaceLinear = svs.createPiecewiseLinearContinuousVectorSpace(grid)
spaceDLinear = svs.createPiecewiseLinearContinuousVectorSpace(grid, segmentD)
spaceNLinear = svs.createPiecewiseLinearContinuousVectorSpace(grid, segmentN)
spaceDLinearClipped = svs.createPiecewiseLinearContinuousVectorSpace(grid, segmentD, True)

spaceDis = svs.createPiecewiseLinearDiscontinuousVectorSpace(grid)

# Integral equations
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Left hand side operator
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
slpOpDD = sln.createNavier3dSingleLayerBoundaryOperator(context, 
                                                    spaceDLinear, spaceDLinear, 
                                                    spaceDLinear, 
                                                    k_p, k_s, mu, "SLP DD")
slpOpND = sln.createNavier3dSingleLayerBoundaryOperator(context, 
                                                    spaceDLinear, spaceNLinear, 
                                                    spaceNLinear, 
                                                    k_p, k_s, mu, "SLP ND")
dlpOpDN = dln.createNavier3dDoubleLayerBoundaryOperator(context, 
                                                    spaceNLinear, spaceDLinear, 
                                                    spaceDLinear, 
                                                    k_p, k_s, "DLP DN")
dlpOpNN = dln.createNavier3dDoubleLayerBoundaryOperator(context, 
                                                    spaceNLinear, spaceNLinear, 
                                                    spaceNLinear, 
                                                    k_p, k_s, "DLP NN")
iOpNN = createIdentityOperator(context, 
                               spaceNLinear, spaceNLinear, 
                               spaceNLinear, 
                               "I NN")

lhsOp = createBlockedBoundaryOperator(context,
                                      [ [  slpOpDD,       -dlpOpDN       ],
                                        [  slpOpND, -(0.5*iOpNN+dlpOpNN) ] ])
# Right hand side
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
slpOpDN = sln.createNavier3dSingleLayerBoundaryOperator(context, 
                                                    spaceNLinear, spaceDLinear, 
                                                    spaceDLinear, 
                                                    k_p, k_s, mu, "SLP DN")
slpOpNN = sln.createNavier3dSingleLayerBoundaryOperator(context, 
                                                    spaceNLinear, spaceNLinear, 
                                                    spaceNLinear, 
                                                    k_p, k_s, mu, "SLP NN")
dlpOpDD = dln.createNavier3dDoubleLayerBoundaryOperator(context, 
                                                    spaceDLinear, spaceDLinear, 
                                                    spaceDLinear, 
                                                    k_p, k_s, "DLP DD")
dlpOpND = dln.createNavier3dDoubleLayerBoundaryOperator(context, 
                                                    spaceDLinear, spaceNLinear, 
                                                    spaceNLinear, 
                                                    k_p, k_s, "DLP ND")
iOpDD = createIdentityOperator(context, 
                               spaceDLinear, spaceDLinear, 
                               spaceDLinear, 
                               "I DD")
# Boundary conditions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def evalDirichletTrace(point):
    x, y, z = point
    return np.array([0., 0., -1.0 + 0j])
dirichletTraceD = createGridFunction(context,spaceDLinear,spaceDLinearClipped,
                                     evalDirichletTrace)

def evalNeumannTrace(point):
    x, y, z = point
    return np.array([0., 0., 0.])
neumannTraceN = createGridFunction(context,spaceNLinear,spaceNLinear,evalNeumannTrace,
                                   surfaceNormalDependent=False)

rhs = [ -slpOpDN*neumannTraceN + (0.5*iOpDD+dlpOpDD)*dirichletTraceD,
        -slpOpNN*neumannTraceN + dlpOpND*dirichletTraceD ]

# Solution of the boundary integral equations
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# solver = createDefaultIterativeSolver(lhsOp)
# solver.initializeSolver(defaultGmresParameterList(1e-8,10000))
solver = createDefaultDirectSolver(lhsOp)

solution = solver.solve(rhs)
print solution.solverMessage()

neumannTraceD = solution.gridFunction(0)
dirichletTraceN = solution.gridFunction(1)

scatterSpaceD = createIdentityOperator(context, spaceDLinear, spaceDis, spaceDis)
scatterSpaceN = createIdentityOperator(context, spaceNLinear, spaceDis, spaceDis)

dirichletTrace = (scatterSpaceD*dirichletTraceD +
                  scatterSpaceN*dirichletTraceN)
neumannTrace = (scatterSpaceD*neumannTraceD +
                scatterSpaceN*neumannTraceN)

dirichletTrace.exportToVtk("vertex_data", "dirichlet_data", "dirichlet_data")
neumannTrace.exportToVtk("vertex_data", "neumann_data", "neumann_data")




