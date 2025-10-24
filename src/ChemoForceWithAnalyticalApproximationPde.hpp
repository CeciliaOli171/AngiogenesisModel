#ifndef CHEMOFORCEWITHANALYTICALAPPROXIMATIONPDE_HPP_
#define CHEMOFORCEWITHANALYTICALAPPROXIMATIONPDE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractLinearPde.hpp"
#include "AbstractBoxDomainPdeModifier.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "VesselTipMutationState.hpp"
#include "VesselSegmentMutationState.hpp"

#include "ChemoForce.hpp"

/**
 * A chemotactic force subclass. Describes the movement of vessel tip due to chemotactic gradient (such as VEGF) and chemoattractant factors' concentrations in angiogenesis. Here, we consider that we are at a steady-state of the PDE system describing the vegf concentration to obtain an analytical approximation formula of vegf concentration.
 */

template<unsigned DIM>
class ChemoForceWithAnalyticalApproximationPde : public ChemoForce<DIM>
{
    friend class TestForcesModel;
    friend class TestAngiogenesisModel;
    friend class TestAngiogenesisModelWithVegfConcentrationAnalyticalApproximationOfPde;

private:
    /* parameters */
    double mChiAnalyticalApproxPde;
    double mHX;
    double mDiffusionCoefficient;
    double mDecayCoefficient;
    double mCreationCoefficient;
    double mConsumptionCoefficient;
    double mSourceValue;
    double mConstantBackground;
    std::vector<c_vector<double, DIM> > mGradientsVegfAnalyticalApproxPde;

    /* serialisation */
    friend class boost::serialization::access;

    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mChiAnalyticalApproxPde;
        archive & mDiffusionCoefficient;
        archive & mDecayCoefficient;
        archive & mCreationCoefficient;
        archive & mConsumptionCoefficient;
        archive & mSourceValue;
        archive & mConstantBackground;
        archive & mGradientsVegfAnalyticalApproxPde;
    }

public:

    /**
     * Constructor.
     *
     * @param chiAnalyticalApproxPde the chemotactic sensitivity in the case of an analytical approximation of vegf concentration
     * @param hx the chemoattractant gradient factor in the x axis
     * @param diffusionCoefficient the vegf diffusion coefficient in the ECM
     * @param decayCoefficient the natural vegf decay
     * @param creationCoefficient the creation of vegf coefficient
     * @param consumptionCoefficient the consumption of vegf by vessel tip coefficient
     * @param sourceValue the value of the source of vegf
     * @param constantBackground the constant vegf concentration value present in the ECM
     */
    ChemoForceWithAnalyticalApproximationPde(double chiAnalyticalApproxPde = 1e-4, double hx=1e-4, double diffusionCoefficient=1.0, double decayCoefficient=1.0, double creationCoefficient=0.1, double consumptionCoefficient=0.01, double sourceValue=0.5, double constantBackground=0.1);

    /**
     * Destructor.
     */
    ~ChemoForceWithAnalyticalApproximationPde();

    /**
     * Overridden GetGradient() method.
     *
     * Get the vegf gradient at a given node.
     *
     * @param rCellPopulation reference to the cell population
     * @param node_index index of the node
     *
     * @return vegf gradient associates to node
     */
    c_vector<double, DIM>& GetGradient(AbstractCellPopulation<DIM>& rCellPopulation, unsigned node_index);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemoForceWithAnalyticalApproximationPde)

#endif /*CHEMOFORCEWITHANALYTICALAPPROXIMATIONPDE_HPP_*/