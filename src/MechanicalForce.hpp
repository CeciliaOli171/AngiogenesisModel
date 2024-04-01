#ifndef MECHANICALFORCE_HPP_
#define MECHANICALFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "SmartPointers.hpp"

template<unsigned DIM>
class MechanicalForce  : public AbstractForce<DIM>
{
friend class TestForces;

private:

    double mSc;
    double GetSpringConstant(); // sensitivity to random fluctuations

    // allow to archive the force model object in a cell-based simulation 
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mSc;
    }

public:

    // constructor 
    MechanicalForce(double Sc=5.56E-6);

    // destructor 
    ~MechanicalForce();

    // calculates the length of a vessel element 
    double LengthVesselSegment(AbstractCellPopulation<DIM, DIM>& rCellPopulation, CellPtr pCell, std::set<unsigned> neighbouring_node_indices);

    // Spring function 
    double SpringFunction(double x, double lc, double Rc);

    // overrides AddForceContribution
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    // overrides OutputForceParameters
    void OutputForceParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MechanicalForce) 

#endif /*MECHANICALFORCE_HPP_*/