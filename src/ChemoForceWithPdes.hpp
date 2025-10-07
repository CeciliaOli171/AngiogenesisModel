#ifndef CHEMOFORCEWITHPDES_HPP_
#define CHEMOFORCEWITHPDES_HPP_

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
 * A chemotactic force subclass. Describes the movement of vessel tip due to chemotactic gradient (such as VEGF) and chemoattractant factors' concentrations in angiogenesis. Here, we consider a reaction-diffusion PDE of vegf concentration.
 */

template<unsigned DIM>
class ChemoForceWithPdes  : public ChemoForce<DIM>
{
    friend class TestForces;
    friend class TestAngiogenesisModel;
    friend class TestAngiogenesisModelWithVegfConcentrationPde;

private:
    /* parameters */
    double mChiPdes;
    boost::shared_ptr<AbstractBoxDomainPdeModifier<DIM> > mpPdeModifier;
    std::vector<c_vector<double, DIM> > mGradientsVegfPdes;

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
        archive & mChiPdes;
        archive & mpPdeModifier;
        archive & mGradientsVegfPdes;
    }

public:
    /**
     * Constructor.
     *
     * @param chiPdes the chemotactic sensitivity in the case of a reaction-diffusion pde for vegf concentration
     * @param hx the chemoattractant gradient factor in the x axis
     * @param hy the chemoattractant gradient factor in the y axis
     * @param hz the chemoattractant gradient factor in the z axis
     * @param pPdeModifier a shared pointer to a PDE modifier
     */
    ChemoForceWithPdes(double chiPdes = 0.1, double hx=1e-4, double hy=1e-4, double hz=1e-4, boost::shared_ptr<AbstractBoxDomainPdeModifier<DIM> > pPdeModifier=boost::shared_ptr<AbstractBoxDomainPdeModifier<DIM> >());

    /**
     * Destructor.
     */
    ~ChemoForceWithPdes();

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
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemoForceWithPdes)

#endif /*CHEMOFORCEWITHPDES_HPP_*/
