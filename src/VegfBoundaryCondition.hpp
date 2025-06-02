#ifndef _VEGFBOUNDARYCONDITION_HPP_
#define _VEGFBOUNDARYCONDITION_HPP_

#include "AbstractBoundaryCondition.hpp"

#include <boost/serialization/base_object.hpp>


template<unsigned SPACE_DIM>
class VegfBoundaryCondition : public AbstractBoundaryCondition<SPACE_DIM>
{
private:

    friend class boost::serialization::access;
    
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractBoundaryCondition<SPACE_DIM> >(*this);
    }

    double mValue;
    double mBoundaryCuboidMax;

public:

    VegfBoundaryCondition(const double value, double boundaryCuboidMax);

    double GetValue(const ChastePoint<SPACE_DIM>& rX) const;
};

#include "SerializationExportWrapper.hpp" // Must be last
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VegfBoundaryCondition)

namespace boost
{
namespace serialization
{

template<class Archive, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const VegfBoundaryCondition<SPACE_DIM> * t, const unsigned int file_version)
{
    const ChastePoint<SPACE_DIM> p;
    const double value = t->GetValue(p);

    ar & value;
}

template<class Archive, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, VegfBoundaryCondition<SPACE_DIM> * t, const unsigned int file_version)
{
    double value;
    double boundaryCuboidMax;
    ar & value;
    ar & boundaryCuboidMax;

    ::new(t)VegfBoundaryCondition<SPACE_DIM>(value, boundaryCuboidMax);
}
}
} 

#endif //_VEGFBOUNDARYCONDITION_HPP_