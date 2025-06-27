#ifndef PERSISTENCEFORCE_HPP_
#define PERSISTENCEFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "TipCellMutationState.hpp"
#include "VesselCellMutationState.hpp"

#include "SmartPointers.hpp"


template<unsigned DIM>
class PersistenceForce  : public AbstractForce<DIM>
{
friend class TestForcesModel;

private:

    double mOmegap;
    double GetPersistenceCoefficient(); // directional persistence coefficient 

    // allow to archive the force model object in a cell-based simulation 
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mOmegap;
    }

public:

    // constructor 
    PersistenceForce(double omegap = 0.4);

    // destructor 
    ~PersistenceForce();

    // overrides AddForceContribution
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    // overrides OutputForceParameters
    void OutputForceParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PersistenceForce) 

#endif /*PERSISTENCEFORCE_HPP_*/