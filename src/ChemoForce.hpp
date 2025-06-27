#ifndef CHEMOFORCE_HPP_
#define CHEMOFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractLinearPde.hpp"
#include "AbstractBoxDomainPdeModifier.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "TipCellMutationState.hpp"
#include "VesselCellMutationState.hpp"

template<unsigned DIM>
class ChemoForce  : public AbstractForce<DIM>
{
friend class TestForcesModel;

private:

    double mChi;
    double mCX;
    double mCY;
    double mCZ;

    double GetChemotacticSensitivity();
    double GetChemotacticGradientCoefficientXAxis();
    double GetChemotacticGradientCoefficientYAxis();
    double GetChemotacticGradientCoefficientZAxis();

    std::vector<c_vector<double, DIM> > mGradients;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
    }

public:

    // constructor
    ChemoForce(double chi = 5.56E-4, double cx = 0.1, double cy = 0.0, double cz = 0.0);

    // destructor
    ~ChemoForce();

    // display gradient at a node position 
    virtual c_vector<double, DIM>& GetGradient(unsigned node_index);

    // display the magnitude of the gradient at a node position 
    virtual double GetMagnitudeGradient(unsigned node_index);

    // calculate the gradient of the vegf concentration 
    virtual void CalculateVegfGradient(AbstractCellPopulation<DIM>& rCellPopulation);

    // override AddForceContribution
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    // override OutputForceParameters
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemoForce)

#endif /*CHEMOFORCE_HPP_*/
