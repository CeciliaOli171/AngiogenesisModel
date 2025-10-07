#ifndef SPROUTINGRULEWITHCONSTANTVEGF_HPP_
#define SPROUTINGRULEWITHCONSTANTVEGF_HPP_

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
 * A sprouting rule subclass for use in centre-based simulations implemented for angiogenesis. The purpose of this class is to return the locations of the two daughters of a dividing cell depending on the type of division. Vessel tips can either divide into either a VS and a VT, lengthening the vessel (asymmetric division); or into two VTs, creating a new sprout (symmetric division). Here, the sprouting probability depend on the vegf concentration around the vessel tip; we consider that the vegf concentration is constant.
 */

// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCentreBasedCellPopulation;
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCentreBasedDivisionRule;

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class SproutingRuleWithConstantVegf : public SproutingRule<ELEMENT_DIM, SPACE_DIM>
{
    friend class TestForcesModel;
    friend class TestAngiogenesisModelWithVegfConcentrationConstant;

private:
    /* parameters */
    double mSproutingProbabilityWithConstantVegf;
    double mMaxSproutingRateConstantVegf;
    double mCMax;
    double mCMin;
    double mPMax;
    double mPMin;
    double mConstantBackground;

    /**
     * @return the sprouting probability in the case of constant vegf
     */
    double GetSproutingProbabilityWithConstantVegf();

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
        // (not in AbstractCentreBasedDivisionRule class
        archive & boost::serialization::base_object<AbstractCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive & mSproutingProbabilityWithConstantVegf;
        archive & mMaxSproutingRateConstantVegf;
        archive & mCMax;
        archive & mCMin;
        archive & mPMax;
        archive & mPMin;
        archive & mConstantBackground;
    }

public:
    /**
     * Constructor.
     * 
     * @param maxSproutingRateConstantVegf the maximum sprouting rate for constant vegf
     * @param thresholdLength the minimum branch length
     * @param constantBackground the constant vegf concentration value present in the ECM
     * @param cMax the maximum vegf concentration
     * @param cMin the minimum vegf concentration
     * @param pMax the maximum probability of sprouting 
     * @param pMin the minimum probability of sprouting 
     */
    SproutingRuleWithConstantVegf(double MaxSproutingRateConstantVegf = 0.08, double thresholdLength=2.0, double constantBackground=0.5, double cMax=1, double cMin=0.3, double pMax=0.98, double pMin=0.4);

    /**
     * Destructor.
     */
    ~SproutingRuleWithConstantVegf();

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
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SproutingRuleWithConstantVegf)

#endif /*SPROUTINGRULEWITHCONSTANTVEGF_HPP_*/