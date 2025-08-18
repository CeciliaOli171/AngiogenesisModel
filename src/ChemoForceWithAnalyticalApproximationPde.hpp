#ifndef CHEMOFORCEWITHANALYTICALAPPROXIMATIONPDE_HPP_
#define CHEMOFORCEWITHANALYTICALAPPROXIMATIONPDE_HPP_

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
class ChemoForceWithAnalyticalApproximationPde : public ChemoForce<DIM>
{
friend class TestForcesModel;

private:

    double mChiAnalyticalApproxPde;

    double mDiffusionCoefficient;
    double mDecayCoefficient;
    double mCreationCoefficient;
    double mConsumptionCoefficient;
    double mSourceValue;
    double mConstantBackground;

    std::vector<c_vector<double, DIM> > mGradientsVegfAnalyticalApproxPde;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
    }

public:

    // constructor
    ChemoForceWithAnalyticalApproximationPde(double chiAnalyticalApproxPde = 1e-4, double hx=1e-4, double diffusionCoefficient=1.0, double decayCoefficient=1.0, double creationCoefficient=0.1, double consumptionCoefficient=0.01, double sourceValue=0.5, double constantBackground=0.1);

    // destructor
    ~ChemoForceWithAnalyticalApproximationPde();

    // display gradient at a node position 
    c_vector<double, DIM>& GetGradient(unsigned node_index);

    // display the magnitude of the gradient at a node position 
    double GetMagnitudeGradient(unsigned node_index);

    // calculate the gradient of the vegf concentration 
    void CalculateVegfGradient(AbstractCellPopulation<DIM>& rCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemoForceWithAnalyticalApproximationPde)

#endif /*CHEMOFORCEWITHANALYTICALAPPROXIMATIONPDE_HPP_*/