#include "VesselSegmentMutationState.hpp"

VesselSegmentMutationState::VesselSegmentMutationState()
    : AbstractCellMutationState(1)
{}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(VesselSegmentMutationState)
