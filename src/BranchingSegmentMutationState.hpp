#ifndef BRANCHINGCELLMUTATIONSTATE_HPP_
#define BRANCHINGCELLMUTATIONSTATE_HPP_

#include "AbstractCellMutationState.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>


class BranchingSegmentMutationState : public AbstractCellMutationState
{
private:
    
    friend class boost::serialization::access;
    
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellMutationState>(*this);
    }

public:
    
    // constructor
    BranchingSegmentMutationState();
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(BranchingSegmentMutationState)

#endif /* BRANCHINGCELLMUTATIONSTATE_HPP_ */
