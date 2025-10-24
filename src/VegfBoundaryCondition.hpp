#ifndef _VEGFBOUNDARYCONDITION_HPP_
#define _VEGFBOUNDARYCONDITION_HPP_

#include "AbstractBoundaryCondition.hpp"

#include <boost/serialization/base_object.hpp>

/**
 * A class for the boundary condition for the reaction-diffusion pde of vegf concentration for our angiogenesis model.
 */

template<unsigned SPACE_DIM>
class VegfBoundaryCondition : public AbstractBoundaryCondition<SPACE_DIM>
{
    friend class TestForcesModel;
    
private:
    /* parameters */
    double mSourceTerm;
    double mConstantBackground;
    double mBoundaryCuboidMax;

    /**
     * @return source term.
     */
    double GetSourceTerm();
    
    /* serialisation */
    friend class boost::serialization::access;

    /**
     * Archive the member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractBoundaryCondition<SPACE_DIM> >(*this);
    }

public:
    /**
     * Constructor.
     *
     * @param sourceterm vegf source term
     * @param boundaryCuboidMax boundary max of cuboid
     */
    VegfBoundaryCondition(const double sourceterm=0.5, double constantBackground=0.1, double boundaryCuboidMax=10);

    /**
     * Overridden GetValue() method.
     * 
     * Computes the boundary condition at a Chaste point.
     * 
     * @param rX The point at which to evaluate the boundary condition
     *
     * @return the value of the boundary condition at a given point.
     */
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
    const double sourceterm = t->GetValue(p);

    ar & sourceterm;
}

template<class Archive, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, VegfBoundaryCondition<SPACE_DIM> * t, const unsigned int file_version)
{
    double sourceterm;
    double boundaryCuboidMax;

    ar & sourceterm;
    ar & boundaryCuboidMax;

    ::new(t)VegfBoundaryCondition<SPACE_DIM>(sourceterm, boundaryCuboidMax);
}
}
} 

#endif //_VEGFBOUNDARYCONDITION_HPP_