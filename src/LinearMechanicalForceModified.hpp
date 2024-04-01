#ifndef LINEARMECHANICALFORCEMODIFIED_HPP_
#define LINEARMECHANICALFORCEMODIFIED_HPP_

#include "GeneralisedLinearSpringForce.hpp"

template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class LinearMechanicalForceModified : public GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>
{
    friend class TestForces;

private:

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractTwoBodyInteractionForce<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    /**
     * Constructor.
     */
    LinearMechanicalForceModified();

    /**
     * Destructor.
     */
    ~LinearMechanicalForceModified();

    /**
     * Return a multiplication factor for the spring constant, which
     * returns a default value of 1.
     *
     * This method may be overridden in subclasses.
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
