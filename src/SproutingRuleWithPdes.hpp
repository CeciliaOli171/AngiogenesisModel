#ifndef SPROUTINGRULEWITHPDES_HPP_
#define SPROUTINGRULEWITHPDES_HPP_

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
class SproutingRuleWithPdes  : public SproutingRule<ELEMENT_DIM, SPACE_DIM>
{
    friend class TestForcesModel;
    friend class TestAngiogenesisModelWithVegfConcentrationPde;

private:

    double mMaxSproutingRatePdes;

    boost::shared_ptr<AbstractBoxDomainPdeModifier<SPACE_DIM> > mpPdeModifier; 

    double mCMax;
    double mCMin;
    double mPMax;
    double mPMin;

    // allow to archive the force model object in a cell-based simulation 
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // not in AbstractCentreBasedDivisionRule.hpp
        archive & boost::serialization::base_object<AbstractCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive & mMaxSproutingRatePdes;

        archive & mpPdeModifier; 

        archive & mCMax;
        archive & mCMin;
        archive & mPMax;
        archive & mPMin;
    }

public:

    // constructor 
    SproutingRuleWithPdes(double MaxSproutingRatePdes = 0.08, double thresholdLength=2.0, boost::shared_ptr<AbstractBoxDomainPdeModifier<SPACE_DIM> > pPdeModifier=boost::shared_ptr<AbstractBoxDomainPdeModifier<SPACE_DIM> >(), double cMax=1, double cMin=0.3, double pMax=0.98, double pMin=0.4);

    // destructor 
    ~SproutingRuleWithPdes();

    // calculates the sprouting probability depending on the vegf concentration around the tip cell 
    double GetSproutingProbability(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pParentCell);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SproutingRuleWithPdes)

#endif /*SPROUTINGRULEWITHPDES_HPP_*/