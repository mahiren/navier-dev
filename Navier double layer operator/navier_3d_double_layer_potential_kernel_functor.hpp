#ifndef navier_3d_double_layer_potential_kernel_functor_hpp
#define navier_3d_double_layer_potential_kernel_functor_hpp

#include "common/scalar_traits.hpp"
#include "fiber/geometrical_data.hpp"

#include <iostream>

template <typename ValueType_>
class Navier3dDoubleLayerPotentialKernelFunctor
{
public:
    typedef ValueType_ ValueType;
    typedef typename Bempp::ScalarTraits<ValueType_>::RealType CoordinateType;

  explicit Navier3dDoubleLayerPotentialKernelFunctor(ValueType waveNumber_p,
                                                     ValueType waveNumber_s) :
    kp(waveNumber_p), 
    ks(waveNumber_s)
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

	ValueType i     = 1i ;
	ValueType ONE   = static_cast<CoordinateType>(1.);
	ValueType TWO   = static_cast<CoordinateType>(2.);
	ValueType THREE = static_cast<CoordinateType>(3.);
	ValueType FOUR  = static_cast<CoordinateType>(4.);
	ValueType NINE  = static_cast<CoordinateType>(9.);

        CoordinateType r_1 = ( trialGeomData.global(0) - testGeomData.global(0) );
        CoordinateType r_2 = ( trialGeomData.global(1) - testGeomData.global(1) );
        CoordinateType r_3 = ( trialGeomData.global(2) - testGeomData.global(2) );
        CoordinateType r2 = r_1*r_1 + r_2*r_2 + r_3*r_3;
        CoordinateType r = sqrt(r2);
	CoordinateType r3 = r*r2;

	CoordinateType n_1 = trialGeomData.normal(0);
	CoordinateType n_2 = trialGeomData.normal(1);
	CoordinateType n_3 = trialGeomData.normal(2);
        CoordinateType r_dot_n = r_1*n_1 + r_2*n_2 + r_3*n_3;

	ValueType ksr   = ks * r;
	ValueType kpr   = kp * r;
	ValueType ksr2  = ksr * ksr;
	ValueType kpr2  = kpr * kpr;    

	ValueType sinc_ks = exp( i * ksr ) / r;
	ValueType sinc_kp = exp( i * kpr ) / r;

	ValueType phi = 0., dphi_dr = 0.;
	ValueType psi = 0., dpsi_dr = 0.;

	phi += ( ONE + i/ksr - ONE/ksr2 ) * sinc_ks;
	phi -= ( i/kpr - ONE/kpr2 ) * sinc_kp;

        dphi_dr = (i*i*exp(i*ks*r))/r2 - (i*i*exp(i*kp*r))/r2 - exp(i*ks*r)/r2 - 
		  (THREE*exp(i*kp*r))/(kpr2*r2) + (THREE*exp(i*ks*r))/(ksr2*r2) + 
        	  (i*ks*exp(i*ks*r))/r + (THREE*i*exp(i*kp*r))/(kp*r3) - 
                  (THREE*i*exp(i*ks*r))/(ks*r3);

	psi += ( ONE + THREE*i/ksr - THREE/ksr2 ) * sinc_ks;
	psi -= ( ONE + THREE*i/kpr - THREE/kpr2 ) * sinc_kp;

	dpsi_dr = exp(i*kp*r)/r2 - exp(i*ks*r)/r2 - (THREE*i*i*exp(i*kp*r))/r2 + 
		  (THREE*i*i*exp(i*ks*r))/r2 - (NINE*exp(i*kp*r))/(kpr2*r2) + 
		  (NINE*exp(i*ks*r))/(ksr2*r2) - (i*kp*exp(i*kp*r))/r + 
		  (i*ks*exp(i*ks*r))/r + (NINE*i*exp(i*kp*r))/(kp*r3) - 
        	  (NINE*i*exp(i*ks*r))/(ks*r3);

        ValueType c = ONE / ( FOUR * M_PI * r );
	ValueType A = dphi_dr - psi/r;
	ValueType B = TWO/r2*( TWO*psi/r - dpsi_dr );
	ValueType C = ks*ks/(kp*kp)*( dphi_dr - dpsi_dr - TWO*psi/r ) - 
                                    TWO*( dphi_dr - dpsi_dr - psi/r );

	result[0](0,0) = c * ( A * ( r_dot_n + n_1 * r_1 ) + B * r_1 * r_1 * r_dot_n + C * r_1 * n_1 );
	result[0](0,1) = c * ( A * (           n_1 * r_2 ) + B * r_1 * r_2 * r_dot_n + C * r_1 * n_2 );
	result[0](0,2) = c * ( A * (           n_1 * r_3 ) + B * r_1 * r_3 * r_dot_n + C * r_1 * n_3 );
	result[0](1,0) = c * ( A * (           n_2 * r_1 ) + B * r_2 * r_1 * r_dot_n + C * r_2 * n_1 );
	result[0](1,1) = c * ( A * ( r_dot_n + n_2 * r_2 ) + B * r_2 * r_2 * r_dot_n + C * r_2 * n_2 );
	result[0](1,2) = c * ( A * (           n_2 * r_3 ) + B * r_2 * r_3 * r_dot_n + C * r_2 * n_3 );
	result[0](2,0) = c * ( A * (           n_3 * r_1 ) + B * r_3 * r_1 * r_dot_n + C * r_3 * n_1 );
	result[0](2,1) = c * ( A * (           n_3 * r_2 ) + B * r_3 * r_2 * r_dot_n + C * r_3 * n_2 );
	result[0](2,2) = c * ( A * ( r_dot_n + n_3 * r_3 ) + B * r_3 * r_3 * r_dot_n + C * r_3 * n_3 );

	// std::cout << kp << std::endl;
	// std::cout << ks << std::endl;
	// std::cout << result[0](0,0) << ", " << result[0](0,1) << ", " << result[0](0,2) << std::endl;
	// std::cout << result[0](1,0) << ", " << result[0](1,1) << ", " << result[0](1,2) << std::endl;
	// std::cout << result[0](2,0) << ", " << result[0](2,1) << ", " << result[0](2,2) << std::endl;
	// std::cout << " " << std::endl;

    }

private:
    ValueType kp;
    ValueType ks;
};

#endif
