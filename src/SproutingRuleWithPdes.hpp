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

/**
 * A sprouting rule subclass for use in centre-based simulations implemented for angiogenesis. The purpose of this class is to return the locations of the two daughters of a dividing cell depending on the type of division. Vessel tips can either divide into either a VS and a VT, lengthening the vessel (asymmetric division); or into two VTs, creating a new sprout (symmetric division). Here, the sprouting probability depend on the vegf concentration around the vessel tip; we consider that the vegf concentration is described by a reaction-diffusion pde.
 */

// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCentreBasedCellPopulation;
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCentreBasedDivisionRule;

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class SproutingRuleWithPdes  : public SproutingRule<ELEMENT_DIM, SPACE_DIM>
{
    friend class TestForcesModel;
    friend class TestAngiogenesisModelWithVegfConcentrationPde;

private:
    /* parameters */
    double mMaxSproutingRatePdes;
    boost::shared_ptr<AbstractBoxDomainPdeModifier<SPACE_DIM> > mpPdeModifier; 
    double mCMax;
    double mCMin;
    double mPMax;
    double mPMin;

    /* serialisation */ 
    friend class boost::serialization::access;

    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
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
    /**
     * Constructor.
     * 
     * @param maxSproutingRateConstantVegf the maximum sprouting rate for a vegf concentration described by a reaction-diffusion pde
     * @param thresholdLength the minimum branch length
     * @param pPdeModifier a shared pointer to a PDE modifier
     * @param cMax the maximum vegf concentration
     * @param cMin the minimum vegf concentration
     * @param pMax the maximum probability of sprouting 
     * @param pMin the minimum probability of sprouting 
     */
    SproutingRuleWithPdes(double MaxSproutingRatePdes = 0.08, double thresholdLength=2.0, boost::shared_ptr<AbstractBoxDomainPdeModifier<SPACE_DIM> > pPdeModifier=boost::shared_ptr<AbstractBoxDomainPdeModifier<SPACE_DIM> >(), double cMax=1, double cMin=0.3, double pMax=0.98, double pMin=0.4);

    /**
     * Destructor.
     */
    ~SproutingRuleWithPdes();

    /**
     * Overridden GetSproutingProbability() method.
     *
     * Computes the probability for a vessel tip to create a new sprout depending on the vegf concentration. 
     * 
     * @param rCellPopulation the centre-based cell population
     * @param pParentCell the cell to divide
     *
     * @return sprouting probability
     */ 
    double GetSproutingProbability(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pParentCell);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SproutingRuleWithPdes)

#endif /*SPROUTINGRULEWITHPDES_HPP_*/