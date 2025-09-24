#include "BranchingSegmentMutationState.hpp"

BranchingSegmentMutationState::BranchingSegmentMutationState()
    : AbstractCellMutationState(2)
{}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(BranchingSegmentMutationState)
