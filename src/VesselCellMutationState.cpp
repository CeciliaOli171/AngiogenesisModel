#include "VesselCellMutationState.hpp"

VesselCellMutationState::VesselCellMutationState()
    : AbstractCellMutationState(1)
{}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(VesselCellMutationState)
