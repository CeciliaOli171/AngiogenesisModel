#ifndef VESSELCELLMUTATIONSTATE_HPP_
#define VESSELCELLMUTATIONSTATE_HPP_

#include "AbstractCellMutationState.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>


class VesselSegmentMutationState : public AbstractCellMutationState
{
    friend class TestForcesModel;

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellMutationState>(*this);
    }

public:
    
    VesselSegmentMutationState();
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(VesselSegmentMutationState)

#endif /* VESSELCELLMUTATIONSTATE_HPP_ */
