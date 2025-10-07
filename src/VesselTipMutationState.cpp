#include "VesselTipMutationState.hpp"

VesselTipMutationState::VesselTipMutationState()
    : AbstractCellMutationState(0)
{}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(VesselTipMutationState)
