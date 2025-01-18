#ifndef CHEMOFORCEWITHPDES_HPP_
#define CHEMOFORCEWITHPDES_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractLinearPde.hpp"
#include "AbstractBoxDomainPdeModifier.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "TipCellMutationState.hpp"
#include "VesselCellMutationState.hpp"

#include "ChemoForce.hpp"

template<unsigned DIM>
class ChemoForceWithPdes  : public ChemoForce<DIM>
{
friend class TestForces;

private:

    double mChi;
    boost::shared_ptr<AbstractBoxDomainPdeModifier<DIM> > mpPdeModifier;

    std::vector<c_vector<double, DIM> > mGradientsWithVegf;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
    }

public:

    // constructor
    ChemoForceWithPdes(double chi = 0.1, boost::shared_ptr<AbstractBoxDomainPdeModifier<DIM> > pPdeModifier=boost::shared_ptr<AbstractBoxDomainPdeModifier<DIM> >());

    // destructor
    ~ChemoForceWithPdes();

    // display gradient at a node position 
    c_vector<double, DIM>& GetGradient(unsigned node_index);

    // calculate the gradient of the vegf concentration 
    void CalculateVegfGradient(AbstractCellPopulation<DIM>& rCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemoForceWithPdes)

#endif /*CHEMOFORCEWITHPDES_HPP_*/
