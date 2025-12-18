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

/**
 * A sprouting rule subclass for use in centre-based simulations implemented for angiogenesis. The purpose of this class is to return the locations of the two daughters of a dividing cell depending on the type of division. Vessel tips can either divide into either a VS and a VT, lengthening the vessel (asymmetric division); or into two VTs, creating a new sprout (symmetric division). Here, the sprouting probability depend on the vegf concentration around the vessel tip; we consider that we are at a steady-state of the PDE system describing the vegf concentration to obtain an analytical approximation formula of vegf concentration.
 */

// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCentreBasedCellPopulation;
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCentreBasedDivisionRule;

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class SproutingRuleWithAnalyticalApproximationPde : public SproutingRule<ELEMENT_DIM, SPACE_DIM>
{
    friend class TestForcesModel;
    friend class TestAngiogenesisModelWithVegfConcentrationAnalyticalApproximationOfPde;

private:
    /* parameters */
    double mMaxSproutingRateAnalyticalApproxPde;
    double mDiffusionCoefficient;
    double mDecayCoefficient;
    double mCreationCoefficient;
    double mConsumptionCoefficient;
    double mCMax;
    double mCMin;
    double mPMax;
    double mPMin;
    double mBoundaryCuboidMax;
    double mSourceValue;
    double mConstantBackground;

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
        // not in (AbstractCentreBasedDivisionRule class)
        archive & boost::serialization::base_object<AbstractCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive & mMaxSproutingRateAnalyticalApproxPde;
        archive & mDiffusionCoefficient;
        archive & mDecayCoefficient;
        archive & mCreationCoefficient;
        archive & mConsumptionCoefficient;
        archive & mCMax;
        archive & mCMin;
        archive & mPMax;
        archive & mPMin;
        archive & mBoundaryCuboidMax;
        archive & mSourceValue;
        archive & mConstantBackground;
    }

public:
    /**
     * Constructor.
     * 
     * @param maxSproutingRateAnalyticalApproxPde the maximum sprouting rate for analytical approximation of vegf concentration
     * @param thresholdLength the minimum branch length
     * @param diffusionCoefficient the vegf diffusion coefficient in the ECM
     * @param decayCoefficient the natural vegf decay
     * @param creationCoefficient the creation of vegf coefficient
     * @param consumptionCoefficient the consumption of vegf by vessel tip coefficient
     * @param sourceValue the value of the source of vegf
     * @param constantBackground the constant vegf concentration value present in the ECM
     * @param cMax the maximum vegf concentration
     * @param cMin the minimum vegf concentration
     * @param pMax the maximum probability of sprouting 
     * @param pMin the minimum probability of sprouting 
     */ 
    SproutingRuleWithAnalyticalApproximationPde(double MaxSproutingRateAnalyticalApproxPde = 0.08, double thresholdLength=2.0, double diffusionCoefficient=1e4, double decayCoefficient=1.0, double creationCoefficient=0.0, double consumptionCoefficient=1.0, double sourceValue=0.5, double constantBackground=0.1, double cMax=1, double cMin=0.3, double pMax=0.98, double pMin=0.4);

    /**
     * Destructor.
     */
    ~SproutingRuleWithAnalyticalApproximationPde();

    /**
     * Computes the vegf concentration at a specific location.
     * 
     * @param rCellPopulation the centre-based cell population
     * @param pParentCell the cell to divide
     *
     * @return the vegf concentration at a specified node
     */
    double GetVegfConcentrationAtNode(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pParentCell);

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
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SproutingRuleWithAnalyticalApproximationPde)

#endif /*SPROUTINGRULEWITHANALYTICALPDEAPPROXIMATION_HPP_*/