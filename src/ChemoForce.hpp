#ifndef CHEMOFORCE_HPP_
#define CHEMOFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractLinearPde.hpp"
#include "AbstractBoxDomainPdeModifier.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "VesselTipMutationState.hpp"
#include "VesselSegmentMutationState.hpp"

/**
 * A chemotactic force class. Describes the movement of vessel tip due to chemotactic gradient (such as VEGF) and chemoattractant factors' concentrations in angiogenesis.
 */

template<unsigned DIM>
class ChemoForce  : public AbstractForce<DIM>
{
friend class TestForcesModel;
friend class TestAngiogenesisModel;
friend class TestAngiogenesisModelWithVegfConcentrationPde;
friend class TestAngiogenesisModelWithVegfConcentrationConstant;
friend class TestAngiogenesisModelWithVegfConcentrationAnalyticalApproximationOfPde;

private:
    /* parameters */
    double mChi;
    double mHX;
    double mHY;
    double mHZ;
    double mCX;
    double mCY;
    double mCZ;
    std::vector<c_vector<double, DIM> > mGradients;

    /**
     * @return the chemotactic sensitivity.
     */
    double GetChemotacticSensitivity();

    /**
     * @return the chemoattractant gradient factor in the x,y,z axes.
     */
    double GetChemoattractantGradientFactorXAxis();
    double GetChemoattractantGradientFactorYAxis();
    double GetChemoattractantGradientFactorZAxis();

    /**
     * @return the chemotactic gradient coefficient in the x,y,z axes.
     */
    double GetChemotacticGradientCoefficientXAxis();
    double GetChemotacticGradientCoefficientYAxis();
    double GetChemotacticGradientCoefficientZAxis();

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
        archive & mChi;
        archive & mHX;
        archive & mHY;
        archive & mHZ;
        archive & mCX;
        archive & mCY;
        archive & mCZ;
    }

public:
    /**
     * Constructor.
     *
     * @param chi the chemotactic sensitivity
     * @param hx the chemoattractant gradient factor in the x axis
     * @param hy the chemoattractant gradient factor in the y axis
     * @param hz the chemoattractant gradient factor in the z axis
     * @param cx the chemotactic gradient coefficient in the x axis
     * @param cy the chemotactic gradient coefficient in the y axis
     * @param cz the chemotactic gradient coefficient in the z axis
     */
    ChemoForce(double chi = 1E-4, double hx = 1E-4, double hy = 0.0, double hz = 0.0, double cx = 0.1, double cy = 0.0, double cz = 0.0);

    /**
     * Destructor.
     */
    ~ChemoForce();

    /**
     * Get the vegf gradient at a given node.
     *
     * @param node_index index of the node
     *
     * @return vegf gradient associates to node
     */
    virtual c_vector<double, DIM>& GetGradient(unsigned node_index);

    /**
     * Get the magnitude of the gradient at a node position.
     * 
     * @param node_index index of the node
     *
     * @return the magnitude of gradient at the node
     */
    virtual double GetMagnitudeGradient(unsigned node_index);

    /**
     * Computes the gradient of the vegf concentration at the nodes.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void CalculateVegfGradient(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the cell population
     *
     * Fh = chi H(nabla c, c, h)
     * (c the VEGF concentration that is time and space dependent, h a fixed gradient of other chemoattractant factors)
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
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemoForce)

#endif /*CHEMOFORCE_HPP_*/
