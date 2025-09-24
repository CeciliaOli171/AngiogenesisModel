#ifndef ANGULARFORCE_HPP_
#define ANGULARFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "VesselSegmentMutationState.hpp"
#include "VesselTipMutationState.hpp"


template<unsigned DIM>
class AngularForce  : public AbstractForce<DIM>
{
friend class TestForcesModel;
friend class TestAngiogenesisModel;
friend class TestAngiogenesisModelWithVegfConcentrationPde;
friend class TestAngiogenesisModelWithVegfConcentrationConstant;
friend class TestAngiogenesisModelWithVegfConcentrationAnalyticalApproximationOfPde;

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

    // constructor (default value according to Perfhal's parameters)
    AngularForce(double omegaa = 1E-5); 

    // destructor 
    ~AngularForce();

    // calculates the nth value of a neighbouring node indices set 
    unsigned GetNthNeighbourIndice(std::set<unsigned> neighbouring_node_indices, int n);

    // selects all the cells in the same branch
    std::set<unsigned> GetSameBranchNeighbours(AbstractCellPopulation<DIM, DIM>& rCellPopulation,NodeBasedCellPopulation<DIM>* p_node_population, CellPtr pCell);

    // calculates the angle between two vectors u and v 
    double GetAngleFromVectors(c_vector<double,DIM> u, c_vector<double,DIM> v);

    // calculates the angle and the two closest nodes to a third one (branching point not considered)
    std::tuple<double, c_vector<double,DIM>, c_vector<double,DIM>> ClosestAngleVesselSegment(AbstractCellPopulation<DIM>& rCellPopulation, CellPtr pCell, std::set<unsigned> neighbouring_node_indices);

    // calculates the smallest angle made by a node with its neighbours 
    std::tuple<double, c_vector<double,DIM>, c_vector<double,DIM>> OptimalAngleVesselElement(AbstractCellPopulation<DIM>& rCellPopulation, CellPtr pCell, std::set<unsigned> neighbouring_node_indices); 

    // overrides AddForceContribution
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    // overrides OutputForceParameters
    void OutputForceParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AngularForce) 

#endif /*ANGULARFORCE_HPP_*/