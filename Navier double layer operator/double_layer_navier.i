%module double_layer_navier
%{
#define SWIG_FILE_WITH_INIT
#include "navier_3d_double_layer_boundary_operator.hpp"
%}

%include "bempp.swg"
%feature("compactdefaultargs") navier3dDoubleLayerBoundaryOperator;
%include "navier_3d_double_layer_boundary_operator.hpp"
%template(createNavier3dDoubleLayerBoundaryOperator) navier3dDoubleLayerBoundaryOperator<double>;

