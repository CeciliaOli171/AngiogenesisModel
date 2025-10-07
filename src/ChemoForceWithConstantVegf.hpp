#ifndef CHEMOFORCEWITHCONSTANTVEGF_HPP_
#define CHEMOFORCEWITHCONSTANTVEGF_HPP_

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
 * A chemotactic force subclass. Describes the movement of vessel tip due to chemotactic gradient (such as VEGF) and chemoattractant factors' concentrations in angiogenesis. Here, VEGF concentration is constant.
 */

template<unsigned DIM>
class ChemoForceWithConstantVegf : public ChemoForce<DIM>
{
    friend class TestForcesModel;
    friend class TestAngiogenesisModel;
    friend class TestAngiogenesisModelWithVegfConcentrationConstant;

private:
    /* parameters */
    double mChiConstantVegf;
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
        archive & mChiConstantVegf;
        archive & mConstantBackground;
        archive & mGradientsVegfAnalyticalApproxPde;
    }

public:
    /**
     * Constructor.
     *
     * @param chiConstantVegf the chemotactic sensitivity in the case of a constant vegf concentration
     * @param hx the chemoattractant gradient factor in the x axis
     * @param constantBackground the constant vegf concentration value present in the ECM
     */
    ChemoForceWithConstantVegf(double chiConstantVegf = 1e-4, double hx=1e-4, double constantBackground = 0.1);

    /**
     * Destructor.
     */
    ~ChemoForceWithConstantVegf();

    /**
     * Overridden GetGradient() method.
     *
     * Get the vegf gradient at a given node.
     *
     * @param node_index index of the node
     *
     * @return vegf gradient associates to node
     */
    c_vector<double, DIM>& GetGradient(unsigned node_index);

    /**
     * Overridden GetMagnitudeGradient() method.
     *
     * Get the magnitude of the gradient at a node position.
     * 
     * @param node_index index of the node
     *
     * @return the magnitude of gradient at the node
     */
    double GetMagnitudeGradient(unsigned node_index);

    /**
     * Overridden CalculateVegfGradient() method.
     *
     * Computes the gradient of the vegf concentration at the nodes.
     *
     * @param rCellPopulation reference to the cell population
     */
    void CalculateVegfGradient(AbstractCellPopulation<DIM>& rCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemoForceWithConstantVegf)

#endif /*CHEMOFORCEWITHCONSTANTVEGF_HPP_*/