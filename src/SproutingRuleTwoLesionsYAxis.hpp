#ifndef SPROUTINGRULETWOLESIONSYAXIS_HPP_
#define SPROUTINGRULETWOLESIONSYAXIS_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/* Header Files */
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "ArchiveOpener.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "CellPropertyCollection.hpp"
#include "CellPropertyRegistry.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "AbstractCentreBasedDivisionRule.hpp"
#include "ClassIsAbstract.hpp"
#include "CellData.hpp"
#include "AbstractLinearPde.hpp"
#include "AbstractBoxDomainPdeModifier.hpp"

#include "SproutingRule.hpp"

// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCentreBasedCellPopulation;
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCentreBasedDivisionRule;

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class SproutingRuleTwoLesionsYAxis : public SproutingRule<ELEMENT_DIM, SPACE_DIM>
{

private:

    double mMaxSproutingRateTwoLesionsYAxis;

    double mDiffusionCoefficient;
    double mDecayCoefficient;
    double mCreationCoefficient;
    double mConsumptionCoefficient;

    double mCMax;
    double mCMin;
    double mPMax;
    double mPMin;

    double mSourceValue1;
    double mSourceValue2;
    double mConstantBackground;

    double mLocationLesion1;
    double mLocationLesion2;

    int mPsproutFunctionTestNb;

    // allow to archive the force model object in a cell-based simulation 
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // not in AbstractCentreBasedDivisionRule.hpp
        archive & boost::serialization::base_object<AbstractCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    // constructor 
    SproutingRuleTwoLesionsYAxis(double MaxSproutingRateTwoLesionsYAxis = 0.08, double thresholdLength=2.0, double diffusionCoefficient=1e4, double decayCoefficient=1.0, double creationCoefficient=0.0, double consumptionCoefficient=1.0, double sourceValue1=0.5, double sourceValue2=0.5, double constantBackground=0.1, double cMax=1, double cMin=0.3, double pMax=0.98, double pMin=0.4, double locationLesion1 = 0.0, double locationLesion2 = 220.0, int PsproutFunctionTestNb = 1);

    // destructor 
    ~SproutingRuleTwoLesionsYAxis();

    // calculates the vegf concentration at a specified node 
    double GetVegfConcentrationAtNode(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pParentCell);

    // calculates the sprouting probability depending on the vegf concentration around the tip cell 
    double GetSproutingProbability(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pParentCell);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SproutingRuleTwoLesionsYAxis)

#endif /*SPROUTINGRULETWOLESIONSYAXIS_HPP_*/