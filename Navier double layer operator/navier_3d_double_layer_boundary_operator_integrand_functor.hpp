// Copyright (C) 2011-2012 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#ifndef navier_3d_double_layer_boundary_operator_integrand_functor_hpp
#define navier_3d_double_layer_boundary_operator_integrand_functor_hpp

#include "../common/common.hpp"

#include "fiber/collection_of_3d_arrays.hpp"
#include "fiber/geometrical_data.hpp"
#include "fiber/conjugate.hpp"

#include <cassert>

namespace Fiber
{

template <typename BasisFunctionType_, typename KernelType_,
          typename ResultType_>
class Navier3dDoubleLayerBoundaryOperatorIntegrandFunctor
{
public:
    typedef BasisFunctionType_ BasisFunctionType;
    typedef KernelType_ KernelType;
    typedef ResultType_ ResultType;
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

    void addGeometricalDependencies(size_t& testGeomDeps, size_t& trialGeomDeps) const {
        // Do nothing
    }

    template <template <typename T> class CollectionOf2dSlicesOfConstNdArrays>
    ResultType evaluate(
            const ConstGeometricalDataSlice<CoordinateType>& testGeomData,
            const ConstGeometricalDataSlice<CoordinateType>& trialGeomData,
            const CollectionOf1dSlicesOfConst3dArrays<BasisFunctionType>& testTransfValues,
            const CollectionOf1dSlicesOfConst3dArrays<BasisFunctionType>& trialTransfValues,
            const CollectionOf2dSlicesOfConstNdArrays<KernelType>& kernelValues) const {
        const int dimWorld = 3;

        // Assert that there are is at least one vector-valued kernel
        assert(kernelValues.size() >= 1);
        assert(kernelValues[0].extent(0) == 3);
        assert(kernelValues[0].extent(1) == 3);

        // Assert that there is at least one test and trial transformation
        // of correct dimensions
        assert(testTransfValues.size() >= 1);
        assert(trialTransfValues.size() >= 1);
        _1dSliceOfConst3dArray<BasisFunctionType> testValues = testTransfValues[0];
        _1dSliceOfConst3dArray<BasisFunctionType> trialValues = trialTransfValues[0];
        assert(testValues.extent(0) == dimWorld);
        assert(trialValues.extent(0) == dimWorld);

	//        uTAv = uT(Av) = uT(A11v1 + A12v2 + A13v3, A21v1 + A22v2 + A23v3, A31v1 + A32v2 + A33v3) = 
	//             u1A11v1 + u1A12v2 + u1A13v3 + 
	//             u2A21v1 + u2A22v2 + u2A23v3 + 
	//             u3A31v1 + u3A32v2 + u3A33v3
	// u = testValues
	// v = trialValues
	// complex conjugate operator (of scalar) is conjugate()

        return testValues(0) * kernelValues[0](0,0) * trialValues(0) +
               testValues(0) * kernelValues[0](0,1) * trialValues(1) +
               testValues(0) * kernelValues[0](0,2) * trialValues(2) +
               testValues(1) * kernelValues[0](1,0) * trialValues(0) +
               testValues(1) * kernelValues[0](1,1) * trialValues(1) +
               testValues(1) * kernelValues[0](1,2) * trialValues(2) +
               testValues(2) * kernelValues[0](2,0) * trialValues(0) +
               testValues(2) * kernelValues[0](2,1) * trialValues(1) +
               testValues(2) * kernelValues[0](2,2) * trialValues(2);

    }
};

} // namespace Fiber

#endif
