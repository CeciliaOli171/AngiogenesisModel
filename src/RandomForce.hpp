#ifndef RANDOMFORCE_HPP_
#define RANDOMFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "SimulationTime.hpp"

/**
 * A random force class. Describes the interactions of the cell with its micro-environment that induce a random movement.
 */

template<unsigned DIM>
class RandomForce  : public AbstractForce<DIM>
{
friend class TestForcesModel;

private:
    /* parameters */
    double mSigma;

    /**
     * @return the sensitivity to random fluctuations.
     */
    double GetRandomSensitivity(); 

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
        archive & mSigma;
    }

public:
    /**
     * Constructor.
     *
     * @param sigma the sensitivity to random fluctuations
     */
    RandomForce(double sigma=0.1);

    /**
     * Destructor.
     */ 
    ~RandomForce();

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the cell population
     *
     * Fr = sigma xi_c
     * (xi random unit vector chosen from a uniform distribution)
     *
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RandomForce)

#endif /*RANDOMFORCE_HPP_*/