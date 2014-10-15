#ifndef navier_3d_double_layer_boundary_operator_hpp
#define navier_3d_double_layer_boundary_operator_hpp

#include "navier_3d_double_layer_potential_kernel_functor.hpp"
#include "navier_3d_double_layer_boundary_operator_integrand_functor.hpp"

#include "assembly/boundary_operator.hpp"
#include "assembly/general_elementary_singular_integral_operator_imp.hpp"
#include "fiber/scalar_function_value_functor.hpp"
#include <boost/make_shared.hpp>

template <typename BasisFunctionType>
Bempp::BoundaryOperator<BasisFunctionType,
    typename Bempp::ScalarTraits<BasisFunctionType>::ComplexType>
navier3dDoubleLayerBoundaryOperator(
        const boost::shared_ptr<const Bempp::Context<BasisFunctionType,
            typename Bempp::ScalarTraits<BasisFunctionType>::ComplexType>
            >& context,
        const boost::shared_ptr<const Bempp::Space<BasisFunctionType> >& domain,
        const boost::shared_ptr<const Bempp::Space<BasisFunctionType> >& range,
        const boost::shared_ptr<const Bempp::Space<BasisFunctionType> >& dualToRange,
        typename Bempp::ScalarTraits<BasisFunctionType>::ComplexType waveNumber_p,
        typename Bempp::ScalarTraits<BasisFunctionType>::ComplexType waveNumber_s,
        const std::string& label = "")
{
    typedef typename Bempp::ScalarTraits<BasisFunctionType>::ComplexType KernelType;
    typedef typename Bempp::ScalarTraits<BasisFunctionType>::ComplexType ResultType;
    typedef typename Bempp::ScalarTraits<BasisFunctionType>::RealType CoordinateType;

    typedef Navier3dDoubleLayerPotentialKernelFunctor<KernelType> KernelFunctor;
    typedef Fiber::ScalarFunctionValueFunctor<CoordinateType>
    TransformationFunctor;
    typedef Fiber::Navier3dDoubleLayerBoundaryOperatorIntegrandFunctor<
    BasisFunctionType, KernelType, ResultType> IntegrandFunctor;

    typedef Bempp::GeneralElementarySingularIntegralOperator<
            BasisFunctionType, KernelType, ResultType> Op;
    return Bempp::BoundaryOperator<BasisFunctionType, ResultType>(
                context, boost::make_shared<Op>(
                    domain, range, dualToRange, label, Bempp::NO_SYMMETRY,
                    KernelFunctor(waveNumber_p,waveNumber_s),
                    TransformationFunctor(),
                    TransformationFunctor(),
                    IntegrandFunctor()));
}

#endif
