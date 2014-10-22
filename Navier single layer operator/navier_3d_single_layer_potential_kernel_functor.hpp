#ifndef navier_3d_single_layer_potential_kernel_functor_hpp
#define navier_3d_single_layer_potential_kernel_functor_hpp

#include "common/scalar_traits.hpp"
#include "fiber/geometrical_data.hpp"

template <typename ValueType_>
class Navier3dSingleLayerPotentialKernelFunctor
{
public:
    typedef ValueType_ ValueType;
    typedef typename Bempp::ScalarTraits<ValueType_>::RealType CoordinateType;

  explicit Navier3dSingleLayerPotentialKernelFunctor(ValueType waveNumber_p,
                                                     ValueType waveNumber_s,
                                                     CoordinateType mu) :
    k_p(ValueType(0., 1.) * waveNumber_p), 
    k_s(ValueType(0., 1.) * waveNumber_s),
    m(mu)
    {}

    int kernelCount() const { return 1; }
    int kernelRowCount(int kernelIndex) const { return 3; }
    int kernelColCount(int kernelIndex) const { return 3; }

    void addGeometricalDependencies(size_t& testGeomDeps,
                                    size_t& trialGeomDeps) const {
        testGeomDeps |= Fiber::GLOBALS;
        trialGeomDeps |= Fiber::GLOBALS | Fiber::NORMALS;
    }

    template <template <typename T> class CollectionOf2dSlicesOfNdArrays>
    void evaluate(
            const Fiber::ConstGeometricalDataSlice<CoordinateType>& testGeomData,
            const Fiber::ConstGeometricalDataSlice<CoordinateType>& trialGeomData,
            CollectionOf2dSlicesOfNdArrays<ValueType>& result) const {
        const int coordCount = 3;
        CoordinateType d_dot_n = 0., d_dot_d = 0.;
        for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex) {
            CoordinateType diff = trialGeomData.global(coordIndex) -
                    testGeomData.global(coordIndex);
            d_dot_d += diff * diff;
            d_dot_n += diff * trialGeomData.normal(coordIndex);
        }
        CoordinateType d = sqrt(d_dot_d);
        CoordinateType grad_d_1 = ( trialGeomData.global(0) - testGeomData.global(0) ) / d;
        CoordinateType grad_d_2 = ( trialGeomData.global(1) - testGeomData.global(1) ) / d;
        CoordinateType grad_d_3 = ( trialGeomData.global(2) - testGeomData.global(2) ) / d;

	ValueType psi = 0., chi = 0.;
	psi += (static_cast<CoordinateType>(1.) - (static_cast<CoordinateType>(1.) + k_s * d) /
	       (-k_s * k_s * d * d)) * exp(-k_s * d) / d;
	psi += (-k_p * k_p) / (k_s * k_s) * (static_cast<CoordinateType>(1.) + k_p * d) /
               (-k_p * k_p * d * d) * exp(-k_p * d) / d;
	chi += static_cast<CoordinateType>(3.) * psi - static_cast<CoordinateType>(2.) * exp(-k_s * d) / d;
        chi -= (-k_p * k_p) / (k_s * k_s) * exp(-k_p * d) / d;

        CoordinateType c = static_cast<CoordinateType>(1.) / 
	  (static_cast<CoordinateType>(4.) * M_PI * static_cast<CoordinateType>(m) );
        result[0](0, 0) = c * (psi - chi * grad_d_1 * grad_d_1);
        result[0](0, 1) = c * (-chi * grad_d_1 * grad_d_2);
        result[0](0, 2) = c * (-chi * grad_d_1 * grad_d_3);
        result[0](1, 0) = c * (-chi * grad_d_2 * grad_d_1);
        result[0](1, 1) = c * (psi - chi * grad_d_2 * grad_d_2);
        result[0](1, 2) = c * (-chi * grad_d_2 * grad_d_3);
        result[0](2, 0) = c * (-chi * grad_d_3 * grad_d_1);
        result[0](2, 1) = c * (-chi * grad_d_3 * grad_d_2);
        result[0](2, 2) = c * (psi - chi * grad_d_3 * grad_d_3);
    }

private:
    ValueType k_p;
    ValueType k_s;
    CoordinateType m;
};

#endif
