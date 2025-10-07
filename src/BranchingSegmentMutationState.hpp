#ifndef BRANCHINGCELLMUTATIONSTATE_HPP_
#define BRANCHINGCELLMUTATIONSTATE_HPP_

#include "AbstractCellMutationState.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * Subclass of AbstractCellMutationState defining a 'branching segment type' mutation state for the angiogenesis model.
 */

class BranchingSegmentMutationState : public AbstractCellMutationState
{
private:
    /* serialisation */
    friend class boost::serialization::access;

    /**
     * Archive the cell mutation state.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellMutationState>(*this);
    }

public:
    /**
     * Constructor.
     */
    BranchingSegmentMutationState();
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(BranchingSegmentMutationState)

#endif /* BRANCHINGCELLMUTATIONSTATE_HPP_ */
