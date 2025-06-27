#ifndef CHEMOFORCEWITHCONSTANTVEGF_HPP_
#define CHEMOFORCEWITHCONSTANTVEGF_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractLinearPde.hpp"
#include "AbstractBoxDomainPdeModifier.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "TipCellMutationState.hpp"
#include "VesselCellMutationState.hpp"

#include "ChemoForce.hpp"

template<unsigned DIM>
class ChemoForceWithConstantVegf : public ChemoForce<DIM>
{
friend class TestForcesModel;

private:

    double mChiConstantVegf;
    double mCX;

    std::vector<c_vector<double, DIM> > mGradientsVegfAnalyticalApproxPde;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
    }

public:

    // constructor
    ChemoForceWithConstantVegf(double chiAnalyticalApproxPde = 1e-4, double cx);

    // destructor
    ~ChemoForceWithConstantVegf();

    // display gradient at a node position 
    c_vector<double, DIM>& GetGradient(unsigned node_index);

    // display the magnitude of the gradient at a node position 
    double GetMagnitudeGradient(unsigned node_index);

    // calculate the gradient of the vegf concentration 
    void CalculateVegfGradient(AbstractCellPopulation<DIM>& rCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemoForceWithConstantVegf)

#endif /*CHEMOFORCEWITHCONSTANTVEGF_HPP_*/