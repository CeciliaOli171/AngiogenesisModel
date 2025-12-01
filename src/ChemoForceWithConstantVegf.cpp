#include "ChemoForceWithConstantVegf.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "VegfEquationPde.hpp"
#include "LinearBasisFunction.hpp"
#include "ReplicatableVector.hpp"

#include "Debug.hpp"

template<unsigned DIM>
ChemoForceWithConstantVegf<DIM>::ChemoForceWithConstantVegf(double chiConstantVegf, double hx)
    : ChemoForce<DIM>(chiConstantVegf, hx, 0.0, 0.0, 0.0, 0.0, 0.0), mChiConstantVegf(chiConstantVegf)
{
    assert(chiConstantVegf>0);
}

template<unsigned DIM>
ChemoForceWithConstantVegf<DIM>::~ChemoForceWithConstantVegf()
{
}

// Explicit instantiation
template class ChemoForceWithConstantVegf<1>;
template class ChemoForceWithConstantVegf<2>;
template class ChemoForceWithConstantVegf<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemoForceWithConstantVegf)