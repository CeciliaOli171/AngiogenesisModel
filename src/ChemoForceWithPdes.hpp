#ifndef CHEMOFORCEWITHPDES_HPP_
#define CHEMOFORCEWITHPDES_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractLinearPde.hpp"
#include "AbstractBoxDomainPdeModifier.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "VesselTipMutationState.hpp"
#include "VesselSegmentMutationState.hpp"

#include "ChemoForce.hpp"


template<unsigned DIM>
class ChemoForceWithPdes  : public ChemoForce<DIM>
{
friend class TestForces;
friend class TestAngiogenesisModel;
friend class TestAngiogenesisModelWithVegfConcentrationPde;

private:

    double mChiPdes;

    boost::shared_ptr<AbstractBoxDomainPdeModifier<DIM> > mpPdeModifier;
    std::vector<c_vector<double, DIM> > mGradientsVegfPdes;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mChiPdes;

        archive & mpPdeModifier;
        archive & mGradientsVegfPdes;
    }

public:

    // constructor
    ChemoForceWithPdes(double chiPdes = 0.1, double hx=1e-4, double hy=1e-4, double hz=1e-4, boost::shared_ptr<AbstractBoxDomainPdeModifier<DIM> > pPdeModifier=boost::shared_ptr<AbstractBoxDomainPdeModifier<DIM> >());

    // destructor
    ~ChemoForceWithPdes();

    // display gradient at a node position 
    c_vector<double, DIM>& GetGradient(unsigned node_index);

    // display the magnitude of the gradient at a node position 
    double GetMagnitudeGradient(unsigned node_index);

    // calculate the gradient of the vegf concentration 
    void CalculateVegfGradient(AbstractCellPopulation<DIM>& rCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemoForceWithPdes)

#endif /*CHEMOFORCEWITHPDES_HPP_*/
