#ifndef CHEMOFORCETWOLESIONSYAXIS_HPP_
#define CHEMOFORCETWOLESIONSYAXIS_HPP_

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
class ChemoForceTwoLesionsYAxis : public ChemoForce<DIM>
{
friend class TestForcesModel;

private:

    double mChiTwoLesionsYAxis;

    double mDiffusionCoefficient;
    double mDecayCoefficient;
    double mCreationCoefficient;
    double mConsumptionCoefficient;
    double mSourceValue1;
    double mSourceValue2;
    double mConstantBackground;

    double mLocationLesion1;
    double mLocationLesion2;

    std::vector<c_vector<double, DIM> > mGradientsVegfTwoLesionsYAxis;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
    }

public:

    // constructor
    ChemoForceTwoLesionsYAxis(double chiTwoLesionsYAxis = 1e-4, double hy = 1e-4, double diffusionCoefficient=1.0, double decayCoefficient=1.0, double creationCoefficient=0.1, double consumptionCoefficient=0.01, double sourceValue1=0.5, double sourceValue2 = 0.5, double constantBackground=0.1, double locationLesion1=0.0, double locationLesion2=220.0);

    // destructor
    ~ChemoForceTwoLesionsYAxis();

    // display gradient at a node position 
    c_vector<double, DIM>& GetGradient(unsigned node_index);

    // display the magnitude of the gradient at a node position 
    double GetMagnitudeGradient(unsigned node_index);

    // calculate the gradient of the vegf concentration 
    void CalculateVegfGradient(AbstractCellPopulation<DIM>& rCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemoForceTwoLesionsYAxis)

#endif /*CHEMOFORCETWOLESIONSYAXIS_HPP_*/