#ifndef TIPCELLMUTATIONSTATE_HPP_
#define TIPCELLMUTATIONSTATE_HPP_

#include "AbstractCellMutationState.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>


class VesselTipMutationState : public AbstractCellMutationState
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellMutationState>(*this);
    }

public:
    
    VesselTipMutationState();
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(VesselTipMutationState)

#endif /* TIPCELLMUTATIONSTATE_HPP_ */
