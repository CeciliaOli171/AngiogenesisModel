#include "VegfBoundaryCondition.hpp"

#include "Debug.hpp"

template<unsigned SPACE_DIM>
VegfBoundaryCondition<SPACE_DIM>::VegfBoundaryCondition(const double sourceterm, double constantBackground, double boundaryCuboidMax)
    : mSourceTerm(sourceterm), mConstantBackground(constantBackground), mBoundaryCuboidMax(boundaryCuboidMax)
{
}

template<unsigned SPACE_DIM>
double VegfBoundaryCondition<SPACE_DIM>::GetSourceTerm(){
    return mSourceTerm;
}

template<unsigned SPACE_DIM>
double VegfBoundaryCondition<SPACE_DIM>::GetValue(const ChastePoint<SPACE_DIM>& rX) const
{
    // for each point, we check if it is located inside the 
    // at this boundary, the vegf concentration is constant and equal to c0=mValue
    if(rX[0] == mBoundaryCuboidMax){
        return -mSourceTerm;
    } else {
        return mConstantBackground;
    }
}

// explicit instantiation 
template class VegfBoundaryCondition<1>;
template class VegfBoundaryCondition<2>;
template class VegfBoundaryCondition<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VegfBoundaryCondition)
