#include "ConnectingCellMutationState.hpp"

ConnectingCellMutationState::ConnectingCellMutationState()
    : AbstractCellMutationState(3)
{}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(ConnectingCellMutationState)