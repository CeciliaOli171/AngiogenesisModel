#include "VegfBoundaryCondition.hpp"

#include "Debug.hpp"

template<unsigned SPACE_DIM>
VegfBoundaryCondition<SPACE_DIM>::VegfBoundaryCondition(const double value, double boundaryCuboidMax)
    : mValue(value), mBoundaryCuboidMax(boundaryCuboidMax)
{
}

template<unsigned SPACE_DIM>
double VegfBoundaryCondition<SPACE_DIM>::GetValue(const ChastePoint<SPACE_DIM>& rX) const
{
    // for each point, we check if it is located at the boundary of the domain (where the endometriotic lesion is) 
    // at this boundary, the vegf concentration is constant and equal to c0=mValue
    if(rX[0] == mBoundaryCuboidMax){
        // periodic boundary condition:
        // obtain time of simulation t
        // mValue*sin(t) ? 
        return mValue;
    } else {
        return 0.0;
    }
}

// explicit instantiation 
template class VegfBoundaryCondition<1>;
template class VegfBoundaryCondition<2>;
template class VegfBoundaryCondition<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VegfBoundaryCondition)
