#ifndef PERSISTENCEFORCE_HPP_
#define PERSISTENCEFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "VesselTipMutationState.hpp"
#include "VesselSegmentMutationState.hpp"

#include "SmartPointers.hpp"

/**
 * A persistence force class. Describes the movement of the vessel tip along the ECM fibres in angiogenesis.
 */

template<unsigned DIM>
class PersistenceForce  : public AbstractForce<DIM>
{
    friend class TestForcesModel;
    friend class TestAngiogenesisModel;
    friend class TestAngiogenesisModelWithVegfConcentrationPde;
    friend class TestAngiogenesisModelWithVegfConcentrationConstant;
    friend class TestAngiogenesisModelWithVegfConcentrationAnalyticalApproximationOfPde;

private:
    /* parameters */
    double mOmegap;

    /**
     * @return the directional persistence coefficient.
     */
    double GetPersistenceCoefficient(); 

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
        archive & mOmegap;
    }

public:
    /**
     * Constructor.
     *
     * @param omegaa the directional persistence coefficient
     */
    PersistenceForce(double omegap = 0.1);

    /**
     * Destructor.
     */  
    ~PersistenceForce();

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the cell population
     *
     * Fp = omegap (x_p(t) - x_p(t-tau))(|x_p(t) - x_p(t-tau)|)
     * (tau directional persistence time)
     * inspired by Perfhal et al. (2017)
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
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PersistenceForce) 

#endif /*PERSISTENCEFORCE_HPP_*/