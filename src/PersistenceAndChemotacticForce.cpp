#include "PersistenceAndChemotacticForce.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "CellData.hpp"
#include "NodeBasedCellPopulation.hpp"

#include "SimulationTime.hpp"

template<unsigned DIM>
PersistenceAndChemotacticForce<DIM>::PersistenceAndChemotacticForce(double cp)
    : AbstractForce<DIM>()
{
    assert(cp > 0);
    mCp = cp;
}

template<unsigned DIM>
PersistenceAndChemotacticForce<DIM>::~PersistenceAndChemotacticForce()
{
}

template<unsigned DIM>
double PersistenceAndChemotacticForce<DIM>::GetPersistenceAndChemotacticCoefficient()
{
    return mCp;
}

// function that returns the closest neighbour of the vessel element k
template<unsigned DIM>
c_vector<double, DIM> PersistenceAndChemotacticForce<DIM>::ClosestNeighbour(AbstractCellPopulation<DIM, DIM>& rCellPopulation, CellPtr pCell, std::set<unsigned> neighbouring_node_indices){
    // we create a new set with the cell's neighbour from the same branch 
    std::set<unsigned> neighbouring_node_indices_same_branch;
    for(std::set<unsigned>::iterator k = neighbouring_node_indices.begin(); k != neighbouring_node_indices.end(); ++k){
        // we collect the cell pointer 
        CellPtr pNeighbourCell = rCellPopulation.GetCellUsingLocationIndex(*k);
        if (pNeighbourCell->GetCellData()->GetItem("BranchLeader") != pCell->GetCellData()->GetItem("BranchLeader")){
            // we add the cell to the new neighbour set 
            neighbouring_node_indices_same_branch.insert(*k);
        }
    }

    // we collect the cell's neighbours from the same branch data (pointers, coordinates)
    CellPtr p_closest_neighbour = rCellPopulation.GetCellUsingLocationIndex(*neighbouring_node_indices_same_branch.begin());
    c_vector<double, DIM> closest_neighbour = rCellPopulation.GetLocationOfCellCentre(p_closest_neighbour);
    double length_closest_neighbour = norm_2(closest_neighbour); 

    // select the element 
    c_vector<double,DIM> xj = rCellPopulation.GetLocationOfCellCentre(pCell);

    // consider the neighbours of the element 
    // we need to check if they are from the same branch or if the closest neighbour is a branching point
    for(std::set<unsigned>::iterator k = neighbouring_node_indices_same_branch.begin();
    k != neighbouring_node_indices_same_branch.end();
    ++k){
        CellPtr p_neighbour_cell_k = rCellPopulation.GetCellUsingLocationIndex(*k);
        c_vector<double, DIM> xk = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_k);

        c_vector<double, DIM> u = rCellPopulation.rGetMesh().GetVectorFromAtoB(xk, xj);
        double length_u = norm_2(u);
    
        if(length_u < length_closest_neighbour){
            length_closest_neighbour = length_u;
            closest_neighbour = xk;
        }
    }

    // depending on the number of neighbours, we calculate the length of the vessel element
    return closest_neighbour;
}

// overrides AddForceContribution()
template<unsigned DIM>
void PersistenceAndChemotacticForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // initialisation 
    c_vector<double, DIM> persistenceandchemotacticforce; 

    NodeBasedCellPopulation<DIM>* p_node_population = dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);
    if (!p_node_population)
    {
        EXCEPTION("NodeBasedCellPopulation only valid for node.");
    }
    
    for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
        node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
        ++node_iter)
    {
        // we initialise the cell data (index and pointer) that we will use later 
        unsigned node_index = (node_iter)->GetIndex();
        CellPtr pCell = rCellPopulation.GetCellUsingLocationIndex(node_index);

        std::set<unsigned> neighbouring_node_indices = p_node_population->GetNeighbouringNodeIndices(node_index);

        if (pCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()) // && dt % 180 == 0 
        {
            // we collect the cell coordinates 
            c_vector<double,DIM> xi = rCellPopulation.GetLocationOfCellCentre(pCell);

            // we look for its closest neighbour 
            c_vector<double, DIM> xj = ClosestNeighbour(rCellPopulation, pCell, neighbouring_node_indices);

            c_vector<double, DIM> u = rCellPopulation.rGetMesh().GetVectorFromAtoB(xi, xj);;
            double magnitude_closestneighbours = norm_2(u);

            if (magnitude_closestneighbours != 0){
                persistenceandchemotacticforce = mCp/magnitude_closestneighbours*u;
            } else {
                persistenceandchemotacticforce = zero_vector<double>(DIM);
            }
            
            rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(persistenceandchemotacticforce);
        }

    }
}

template<unsigned DIM>
void PersistenceAndChemotacticForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<Directional Persistence Coefficient>" << mCp << "</Cp>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class PersistenceAndChemotacticForce<1>;
template class PersistenceAndChemotacticForce<2>;
template class PersistenceAndChemotacticForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PersistenceAndChemotacticForce)