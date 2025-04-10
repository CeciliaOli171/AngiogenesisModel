#ifndef SPROUTINGRULEWITHANALYTICALAPPROXIMATIONPDE_HPP_
#define SPROUTINGRULEWITHANALYTICALAPPROXIMATIONPDE_HPP_

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
class SproutingRuleWithAnalyticalApproximationPde : public SproutingRule<ELEMENT_DIM, SPACE_DIM>
{

private:

    double mMaxSproutingRateAnalyticalApproxPde;

    double mDiffusionCoefficient;
    double mDecayCoefficient;
    double mCreationCoefficient;
    double mConsumptionCoefficient;

    double mBoundaryCuboidMax;
    double mSourceValue;
    double mConstantBackground;

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
    SproutingRuleWithAnalyticalApproximationPde(double MaxSproutingRateAnalyticalApproxPde = 0.1, double diffusionCoefficient=1.0, double decayCoefficient=1.0, double creationCoefficient=0.1, double consumptionCoefficient=0.01, double boundaryCuboidMax=20.0, double sourceValue=0.1, double constantBackground=0.1, int PsproutFunctionTestNb = 0);

    // destructor 
    ~SproutingRuleWithAnalyticalApproximationPde();

    // calculates the vegf concentration at a specified node 
    double GetVegfConcentrationAtNode(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pParentCell);

    // calculates the sprouting probability depending on the vegf concentration around the tip cell 
    double GetSproutingProbability(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pParentCell);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SproutingRuleWithAnalyticalApproximationPde)

#endif /*SPROUTINGRULEWITHANALYTICALPDEAPPROXIMATION_HPP_*/