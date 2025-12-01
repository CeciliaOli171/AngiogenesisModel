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
    double mHX;
    double mHY;
    double mHZ;
    boost::shared_ptr<AbstractPdeModifier<DIM> > mpPdeModifier;
    std::vector<c_vector<double, DIM> > mGradientsOnElements;
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
    ChemoForceWithPdes(double chiPdes = 0.1, double hx=1e-4, double hy=1e-4, double hz=1e-4, boost::shared_ptr<AbstractPdeModifier<DIM> > pPdeModifier=boost::shared_ptr<AbstractPdeModifier<DIM> >());

    /**
     * Destructor.
     */
    ~ChemoForceWithPdes();

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
    c_vector<double, DIM> GetGradient(AbstractCellPopulation<DIM>& rCellPopulation, unsigned node_index);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemoForceWithPdes)

#endif /*CHEMOFORCEWITHPDES_HPP_*/
