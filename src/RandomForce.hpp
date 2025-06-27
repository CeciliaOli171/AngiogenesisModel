#ifndef RANDOMFORCE_HPP_
#define RANDOMFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "SimulationTime.hpp"

template<unsigned DIM>
class RandomForce  : public AbstractForce<DIM>
{
friend class TestForcesModel;

private:

    double mSigma;
    double GetRandomSensitivity(); // sensitivity to random fluctuations

    // allow to archive the force model object in a cell-based simulation 
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mSigma;
    }

public:

    // constructor 
    RandomForce(double sigma=0.4);

    // destructor 
    ~RandomForce();

    // overrides AddForceContribution
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    // overrides OutputForceParameters
    void OutputForceParameters(out_stream& rParamsFile);

};

// #include "SerializationExportWrapper.hpp"
// EXPORT_TEMPLATE_CLASS_SAME_DIMS(RandomForce)

#endif /*RANDOMFORCE_HPP_*/