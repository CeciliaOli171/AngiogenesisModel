#ifndef PERSISTENCEANDCHEMOTACTICFORCE_HPP_
#define PERSISTENCEANDCHEMOTACTICFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "SmartPointers.hpp"


template<unsigned DIM>
class PersistenceAndChemotacticForce  : public AbstractForce<DIM>
{
friend class TestForces;

private:

    double mCp;
    double GetPersistenceAndChemotacticCoefficient(); // directional persistence coefficient 

    // allow to archive the force model object in a cell-based simulation 
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mCp;
    }

public:

    // constructor 
    PersistenceAndChemotacticForce(double cp = 0.4);

    // destructor 
    ~PersistenceAndChemotacticForce();

    // calculates the closest neighbour to the tip cell 
    c_vector<double, DIM> ClosestNeighbour(AbstractCellPopulation<DIM, DIM>& rCellPopulation, CellPtr pCell, std::set<unsigned> neighbouring_node_indices);

    // overrides AddForceContribution
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    // overrides OutputForceParameters
    void OutputForceParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PersistenceAndChemotacticForce) 

#endif /*PERSISTENCEFORCE_HPP_*/