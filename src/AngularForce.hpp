#ifndef ANGULARFORCE_HPP_
#define ANGULARFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "StemCellProliferativeType.hpp"
#include "AbstractCellProliferativeType.hpp"


template<unsigned DIM>
class AngularForce  : public AbstractForce<DIM>
{
friend class TestForces;

private:

    double mOmegaa;
    double GetAngularPersistence(); // angular persistence spring constant 

    // allow to archive the force model object in a cell-based simulation 
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mOmegaa;
    }

public:

    // constructor 
    AngularForce(double omegaa = 5.56E1); // 5.56E-5 for Perfhal's paper

    // destructor 
    ~AngularForce();

    // smallest branching angle between vessel elements
    std::tuple<double, c_vector<double,DIM>, c_vector<double,DIM>> GetAngleVesselElement(CellPtr cell_ptr, std::set<unsigned> neighbouring_node_indices, AbstractCellPopulation<DIM>& rCellPopulation); 

    // overrides AddForceContribution
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    // overrides OutputForceParameters
    void OutputForceParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AngularForce) 
//CHASTE_CLASS_EXPORT(AngularForce)

#endif /*ANGULARFORCE_HPP_*/