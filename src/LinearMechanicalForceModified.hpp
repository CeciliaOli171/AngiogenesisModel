#ifndef LINEARMECHANICALFORCEMODIFIED_HPP_
#define LINEARMECHANICALFORCEMODIFIED_HPP_

#include "GeneralisedLinearSpringForce.hpp"

/**
 * A linear mechanical force class. Describes the interactions between the cells in a vascular tree in the case of angiogenesis.
 */

template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class LinearMechanicalForceModified : public GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>
{
    friend class TestForcesModel;
    friend class TestAngiogenesisModel;
    friend class TestAngiogenesisModelWithVegfConcentrationPde;
    friend class TestAngiogenesisModelWithVegfConcentrationConstant;
    friend class TestAngiogenesisModelWithVegfConcentrationAnalyticalApproximationOfPde;

private:

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
        archive & boost::serialization::base_object<GeneralisedLinearSpringForce<ELEMENT_DIM,SPACE_DIM> >(*this);
    }

public:
    /**
     * Constructor.
     *
     * @param omegaa the angular persistence spring constant
     */
    LinearMechanicalForceModified();

    /**
     * Destructor.
     */ 
    ~LinearMechanicalForceModified();

    /**
     * Return a multiplication factor for the spring constant, which returns a default value of 1.
     * Applied to angiogenesis: return 0 if two cells are not in the same branch.
     *
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the cell population
     * @param isCloserThanRestLength whether the neighbouring nodes lie closer than the rest length of their connecting spring
     *
     * @return the multiplication factor.
     */
    double VariableSpringConstantMultiplicationFactor(unsigned nodeAGlobalIndex,unsigned nodeBGlobalIndex, AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, bool isCloserThanRestLength);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(LinearMechanicalForceModified)

#endif /*LINEARMECHANICALFORCEMODIFIED_HPP_*/
