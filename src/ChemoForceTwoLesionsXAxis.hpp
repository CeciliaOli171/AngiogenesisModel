#ifndef CHEMOFORCETWOLESIONSXAXIS_HPP_
#define CHEMOFORCETWOLESIONSXAXIS_HPP_

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
class ChemoForceTwoLesionsXAxis : public ChemoForce<DIM>
{
friend class TestForcesModel;

private:

    double mChiTwoLesionsXAxis;

    double mDiffusionCoefficient;
    double mDecayCoefficient;
    double mCreationCoefficient;
    double mConsumptionCoefficient;

    double mSourceValue1;
    double mSourceValue2;
    double mConstantBackground;

    double mLocationLesion1;
    double mLocationLesion2;

    std::vector<c_vector<double, DIM> > mGradientsVegfTwoLesionsXAxis;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
    }

public:

    // constructor
    ChemoForceTwoLesionsXAxis(double chiTwoLesionsXAxis = 1e-4, double hx = 1e-4, double diffusionCoefficient=1.0, double decayCoefficient=1.0, double creationCoefficient=0.1, double consumptionCoefficient=0.01, double sourceValue1=0.5, double sourceValue2=0.5, double constantBackground=0.1, double locationLesion1=55.0, double locationLesion2=165.0);

    // destructor
    ~ChemoForceTwoLesionsXAxis();

    // display gradient at a node position 
    c_vector<double, DIM>& GetGradient(unsigned node_index);

    // display the magnitude of the gradient at a node position 
    double GetMagnitudeGradient(unsigned node_index);

    // calculate the gradient of the vegf concentration 
    void CalculateVegfGradient(AbstractCellPopulation<DIM>& rCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemoForceTwoLesionsXAxis)

#endif /*CHEMOFORCETWOLESIONSXAXIS_HPP_*/