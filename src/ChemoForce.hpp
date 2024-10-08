#ifndef CHEMOFORCE_HPP_
#define CHEMOFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "TipCellMutationState.hpp"
#include "VesselCellMutationState.hpp"

template<unsigned DIM>
class ChemoForce  : public AbstractForce<DIM>
{
friend class TestForces;

private:

    double mChi;
    double mCX;
    double mCY;
    double GetChemotacticSensitivity();
    double GetChemotacticGradientCoefficientXAxis();
    double GetChemotacticGradientCoefficientYAxis();

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
    }

public:

    // constructor
    ChemoForce(double chi = 0.1, double cx = 5.56E-4, double cy = 5.56E-4);

    // destructor
    ~ChemoForce();

    // override AddForceContribution
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    // override OutputForceParameters
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemoForce)

#endif /*CHEMOFORCE_HPP_*/
