#include "BranchingSegmentMutationState.hpp"

BranchingSegmentMutationState::BranchingSegmentMutationState()
    : AbstractCellMutationState(2)
{}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(BranchingSegmentMutationState)
