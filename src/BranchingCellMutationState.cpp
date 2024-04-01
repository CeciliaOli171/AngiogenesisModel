#include "BranchingCellMutationState.hpp"

BranchingCellMutationState::BranchingCellMutationState()
    : AbstractCellMutationState(2)
{}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(BranchingCellMutationState)
