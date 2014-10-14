%module single_layer_navier
%{
#define SWIG_FILE_WITH_INIT
#include "navier_3d_single_layer_boundary_operator.hpp"
%}

%include "bempp.swg"
%feature("compactdefaultargs") navier3dSingleLayerBoundaryOperator;
%include "navier_3d_single_layer_boundary_operator.hpp"
%template(createNavier3dSingleLayerBoundaryOperator) navier3dSingleLayerBoundaryOperator<double>;

