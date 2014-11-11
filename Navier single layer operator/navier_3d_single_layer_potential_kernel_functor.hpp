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
    kp(waveNumber_p), 
    ks(waveNumber_s),
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

	ValueType i     = 1i;
	ValueType ONE   = static_cast<CoordinateType>(1.);
	ValueType THREE = static_cast<CoordinateType>(3.);
	ValueType FOUR  = static_cast<CoordinateType>(4.);

        CoordinateType r_1 = ( trialGeomData.global(0) - testGeomData.global(0) );
        CoordinateType r_2 = ( trialGeomData.global(1) - testGeomData.global(1) );
        CoordinateType r_3 = ( trialGeomData.global(2) - testGeomData.global(2) );
        CoordinateType r = sqrt(r_1*r_1 + r_2*r_2 + r_3*r_3);

	ValueType ksr   = ks * r;
	ValueType kpr   = kp * r;
	ValueType ksr2  = ksr * ksr;
	ValueType kpr2  = kpr * kpr;    

	ValueType sinc_ks = exp( i * ksr ) / r;
	ValueType sinc_kp = exp( i * kpr ) / r;

	ValueType phi = 0.;
	ValueType psi = 0.;

	phi += ( ONE + i/ksr - ONE/ksr2 ) * sinc_ks;
	phi -= ( i/kpr - ONE/kpr2 ) * sinc_kp;

	psi += ( ONE + THREE*i/ksr - THREE/ksr2 ) * sinc_ks;
	psi -= ( ONE + THREE*i/kpr - THREE/kpr2 ) * sinc_kp;

        ValueType c = ONE / ( FOUR * M_PI * static_cast<CoordinateType>(m) );

        result[0](0, 0) = c * ( phi - psi * r_1 * r_1 );
        result[0](0, 1) = c * (     - psi * r_1 * r_2 );
        result[0](0, 2) = c * (     - psi * r_1 * r_3 );
        result[0](1, 0) = c * (     - psi * r_2 * r_1 );
        result[0](1, 1) = c * ( phi - psi * r_2 * r_2 );
        result[0](1, 2) = c * (     - psi * r_2 * r_3 );
        result[0](2, 0) = c * (     - psi * r_3 * r_1 );
        result[0](2, 1) = c * (     - psi * r_3 * r_2 );
        result[0](2, 2) = c * ( phi - psi * r_3 * r_3 );

	// std::cout << kp << std::endl;
	// std::cout << ks << std::endl;
	// std::cout << m << std::endl;
	// std::cout << sinc_ks << "," << sinc_kp << std::endl;
	// std::cout << result[0](0,0) << ", " << result[0](0,1) << ", " << result[0](0,2) << std::endl;
	// std::cout << result[0](1,0) << ", " << result[0](1,1) << ", " << result[0](1,2) << std::endl;
	// std::cout << result[0](2,0) << ", " << result[0](2,1) << ", " << result[0](2,2) << std::endl;
	// std::cout << " " << std::endl;

    }

private:
    ValueType kp;
    ValueType ks;
    CoordinateType m;
};

#endif
