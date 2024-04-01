#include "TipCellMutationState.hpp"

TipCellMutationState::TipCellMutationState()
    : AbstractCellMutationState(0)
{}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(TipCellMutationState)
