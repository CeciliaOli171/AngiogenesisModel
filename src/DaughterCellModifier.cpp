#include "DaughterCellModifier.hpp"
#include "MeshBasedCellPopulation.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "CellData.hpp"

#include "TipCellMutationState.hpp"
#include "VesselCellMutationState.hpp"
#include "BranchingCellMutationState.hpp"

#include "RandomNumberGenerator.hpp"
#include <algorithm>

#include "Debug.hpp" 
#include "SimulationTime.hpp"

template<unsigned DIM>
DaughterCellModifier<DIM>::DaughterCellModifier(double anastomosisLength)
    : AbstractCellBasedSimulationModifier<DIM>(), mHighestBranch(0), mHighestLoop(0), mAnastomosisLength(anastomosisLength)
{
}

template<unsigned DIM>
DaughterCellModifier<DIM>::~DaughterCellModifier()
{
}

template<unsigned DIM>
void DaughterCellModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
std::set<unsigned> DaughterCellModifier<DIM>::GetAnastomosisNeighbours(AbstractCellPopulation<DIM, DIM>& rCellPopulation,NodeBasedCellPopulation<DIM>* p_node_population, CellPtr pParentCell){  
    std::set<unsigned> neighbours_set = p_node_population->GetNodesWithinNeighbourhoodRadius(pParentCell->GetCellId(),mAnastomosisLength);
    std::set<unsigned> neighbours_anastomosis_set;

    for(std::set<unsigned>::iterator k = neighbours_set.begin(); k != neighbours_set.end(); ++k){
        // we collect the cell pointer 
        CellPtr pNeighbourCell = rCellPopulation.GetCellUsingLocationIndex(*k);
        if (pNeighbourCell->GetCellData()->GetItem("BranchNumber") != pParentCell->GetCellData()->GetItem("BranchNumber")){
            // we add the cell to the new neighbour set 
            neighbours_anastomosis_set.insert(*k);
        } 
    }
    return neighbours_anastomosis_set;
}

template<unsigned DIM>
bool DaughterCellModifier<DIM>::IsBranchingCellNextToCell(AbstractCellPopulation<DIM, DIM>& rCellPopulation,NodeBasedCellPopulation<DIM>* p_node_population, CellPtr pParentCell)
{  
    // settings: minimum distance between cell and a branching cell 
    double l_threshold = 2.0;

    // coordinates of parent cell 
    unsigned index_parentcell = rCellPopulation.GetLocationIndexUsingCell(pParentCell);
    c_vector<double, DIM> x_parentcell = rCellPopulation.GetLocationOfCellCentre(pParentCell);

    // get nodes within a given radius 
    std::vector<unsigned>& neighbours_vector = rCellPopulation.GetNode(index_parentcell)->rGetNeighbours(); 
    std::vector<unsigned> neighbours_threshold_vector;

    // we go through all the cells in the neighbours vector
    for(std::vector<unsigned>::iterator i = neighbours_vector.begin(); i != neighbours_vector.end(); ++i)
    {
        // collect the coordinates of the neighbour node 
        CellPtr pNeighbourCell = rCellPopulation.GetCellUsingLocationIndex(*i);
        c_vector<double, DIM> x_neighbourcell = rCellPopulation.GetLocationOfCellCentre(pNeighbourCell);

        // calculate the distance between the nodes 
        c_vector<double, DIM> x_parentneighbourcells = rCellPopulation.rGetMesh().GetVectorFromAtoB(x_parentcell, x_neighbourcell);
        double d_parentneighbourcells = norm_2(x_parentneighbourcells);

        if(d_parentneighbourcells < l_threshold){
            if (pNeighbourCell->GetMutationState()->IsType<BranchingCellMutationState>()){
                return true;
            } 
        }
    }
    return false;
}

// function that returns the closest neighbour of the vessel element k
template<unsigned DIM>
std::pair<c_vector<double, DIM>, unsigned> DaughterCellModifier<DIM>::ClosestNeighbour(AbstractCellPopulation<DIM, DIM>& rCellPopulation, CellPtr pCell, std::set<unsigned> neighbouring_node_indices){
    // we collect the cell's neighbours data (pointers, coordinates)
    CellPtr p_closest_neighbour = rCellPopulation.GetCellUsingLocationIndex(*neighbouring_node_indices.begin());
    c_vector<double, DIM> closest_neighbour = rCellPopulation.GetLocationOfCellCentre(p_closest_neighbour); 
    unsigned indice_closest_neighbour = *neighbouring_node_indices.begin();

    // select the element 
    c_vector<double,DIM> xj = rCellPopulation.GetLocationOfCellCentre(pCell);

    c_vector<double, DIM> u = rCellPopulation.rGetMesh().GetVectorFromAtoB(closest_neighbour, xj);
    double length_closest_neighbour = norm_2(u);

    // consider the neighbours of the element 
    for(std::set<unsigned>::iterator k = neighbouring_node_indices.begin();
    k != neighbouring_node_indices.end();
    ++k){
        CellPtr p_neighbour_cell_k = rCellPopulation.GetCellUsingLocationIndex(*k);
        c_vector<double, DIM> xk = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_k);

        c_vector<double, DIM> u = rCellPopulation.rGetMesh().GetVectorFromAtoB(xk, xj);
        double length_u = norm_2(u);

        // check if it is a branching point or in the same branch 
        if(p_neighbour_cell_k->GetCellData()->GetItem("BranchNumber") == pCell->GetCellData()->GetItem("BranchNumber")){
            if(length_u < length_closest_neighbour){
            length_closest_neighbour = length_u;
            closest_neighbour = xk;
            indice_closest_neighbour = *k;
            }
        } else if(p_neighbour_cell_k->GetMutationState()->IsType<BranchingCellMutationState>() && pCell->GetCellData()->GetItem("BranchingCell") == *k){
            if(length_u < length_closest_neighbour){
            length_closest_neighbour = length_u;
            closest_neighbour = xk;
            indice_closest_neighbour = *k;
            }
        }
    }

    // depending on the number of neighbours, we calculate the length of the vessel element
    return std::make_pair(closest_neighbour, indice_closest_neighbour);
}

template<unsigned DIM>
void DaughterCellModifier<DIM>::CalculateAnastomosisVector( AbstractCellPopulation<DIM, DIM>& rCellPopulation, NodeBasedCellPopulation<DIM>* p_node_population, CellPtr pParentCell){
    // we check if we are in an anastomosis case i.e. two cells with a distance smaller than a cut-off length lanastomosis
    // we select all the neighbours within a neighbourhood radius equal to the cut-off length lanastomosis
    // we set lanastomosis = lc/2 for now

    // consider a parent cell 
    c_vector<double, DIM> parent_position = rCellPopulation.GetLocationOfCellCentre(pParentCell);

    std::set<unsigned> neighbours_anastomosis_set = GetAnastomosisNeighbours(rCellPopulation, p_node_population, pParentCell);
    std::vector<unsigned> neighbours_anastomosis_indices = std::vector<unsigned>(neighbours_anastomosis_set.begin(), neighbours_anastomosis_set.end());

    // select the closest neighbour among the overlapping neighbours 
    unsigned closest_indice;
    if(neighbours_anastomosis_indices.size() >= 1){

        if(neighbours_anastomosis_indices.size() > 1){
            c_vector<double, DIM> closest_neighbour;
            std::tie(closest_neighbour, closest_indice) = ClosestNeighbour(rCellPopulation, pParentCell, neighbours_anastomosis_set);
        } else {
            closest_indice = *neighbours_anastomosis_indices.begin(); // because it is supposed to have only one neighbour in the same branch 
        }
        
        // we collect the closest neighbour data (cell pointer, coordinates, indice)
        CellPtr pClosestNeighbour = rCellPopulation.GetCellUsingLocationIndex(closest_indice);
        c_vector<double, DIM> closest_neighbour_position = rCellPopulation.GetLocationOfCellCentre(pClosestNeighbour);

        // we need to check if the neighbour cell is a vessel segment or a tip cell 
        if (pClosestNeighbour->GetMutationState()->IsType<TipCellMutationState>()){
            // if it is a tip cell, then they merge and become two vessel segment 
            MAKE_PTR(VesselCellMutationState, p_vessel_state);
            pParentCell->SetMutationState(p_vessel_state); // parent cell becomes a vessel segment
            pClosestNeighbour->SetMutationState(p_vessel_state); // closest neighbour becomes a vessel segment 
            MAKE_PTR(DifferentiatedCellProliferativeType, p_vessel_type);
            pParentCell->SetCellProliferativeType(p_vessel_type);
            pParentCell->InitialiseCellCycleModel();
            pClosestNeighbour->SetCellProliferativeType(p_vessel_type);
            pClosestNeighbour->InitialiseCellCycleModel();

            ++mHighestLoop;
            pParentCell->GetCellData()->SetItem("LoopNumber", mHighestLoop);
            pClosestNeighbour->GetCellData()->SetItem("LoopNumber", mHighestLoop);

            // TRACE("Anastomosis between two tip cells");
        } else if (pClosestNeighbour->GetMutationState()->IsType<VesselCellMutationState>()) {
            // if it is a vessel segment, the parent cell becomes a vessel segment and the neighbour cell becomes a branching cell 
            MAKE_PTR(VesselCellMutationState, p_vessel_state);
            pParentCell->SetMutationState(p_vessel_state); // parent cell becomes a vessel segment
            MAKE_PTR(BranchingCellMutationState, p_branching_state);
            pClosestNeighbour->SetMutationState(p_branching_state); // closest neighbour becomes a branching point 

            MAKE_PTR(DifferentiatedCellProliferativeType, p_vessel_type);
            pParentCell->SetCellProliferativeType(p_vessel_type);
            pParentCell->InitialiseCellCycleModel();
            pClosestNeighbour->SetCellProliferativeType(p_vessel_type);
            pClosestNeighbour->InitialiseCellCycleModel();

            pClosestNeighbour->GetCellData()->SetItem("BranchingPoint", 1.0); // the neighour cell becomes a branching cell 
            pParentCell->GetCellData()->SetItem("BranchingLeader", closest_indice);

            ++mHighestLoop;
            pParentCell->GetCellData()->SetItem("LoopNumber", mHighestLoop);
            pClosestNeighbour->GetCellData()->SetItem("LoopNumber", mHighestLoop);

            // TRACE("Anastomosis between a tip cell and a vessel segment");
        } else {
            // if it is a branching cell, the parent cell becomes a vessel segment and the neighbour cell does not change 
            MAKE_PTR(VesselCellMutationState, p_vessel_state);
            pParentCell->SetMutationState(p_vessel_state); // closest neighbour becomes a branching point 

            MAKE_PTR(DifferentiatedCellProliferativeType, p_vessel_type);
            pParentCell->SetCellProliferativeType(p_vessel_type);
            pParentCell->InitialiseCellCycleModel();

            pParentCell->GetCellData()->SetItem("BranchingLeader", closest_indice);

            ++mHighestLoop;
            pParentCell->GetCellData()->SetItem("LoopNumber", mHighestLoop);
            pClosestNeighbour->GetCellData()->SetItem("LoopNumber", mHighestLoop);

            // TRACE("Anastomosis between a tip cell and a branching cell");
        }

        // search for the cells with the same branch number as the parent or neihbour cells 
        // if anastomosis between tip/tip cells or tip/branching cells then the those cells are in the same loop
        // if anastomosis between tip/vessel segment then look for their coordinates and check that they have a smaller norm? 

        double tortuosity_neighbour = pClosestNeighbour->GetCellData()->GetItem("TortuosityParent");
        double tortuosity_parent = tortuosity_neighbour + norm_2(parent_position-closest_neighbour_position);
        pParentCell->GetCellData()->SetItem("TortuosityParent", tortuosity_parent);

        TRACE("ANASTOMOSIS");
    } 
}

template<unsigned DIM>
void DaughterCellModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    // Set up cell data on the cell population : initialisation for the division 
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
    cell_iter != rCellPopulation.End();
    ++cell_iter)
    {
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(node_index);
        c_vector<double, DIM> node_coordinates = rCellPopulation.GetLocationOfCellCentre(p_cell);

        double tortuosity = norm_2(node_coordinates);

        double t = SimulationTime::Instance()->GetTime(); 

        cell_iter->GetCellData()->SetItem("node_index", node_index);
        cell_iter->GetCellData()->SetItem("daughter_type", 0.0); // 0.0 = no daughter cell, 1.0 = daughter cell
        cell_iter->GetCellData()->SetItem("BranchNumber", 0.0);
        cell_iter->GetCellData()->SetItem("BranchingPoint", 0.0); // 0.0 = no, 1.0 = yes : necessary to modify the mutation state after sprouting in this modifier 
        cell_iter->GetCellData()->SetItem("DivisionNumber", 0.0);
        cell_iter->GetCellData()->SetItem("BranchingLeader", 0.0); // indice of the branching cell leading the branch of the cell considered 
        cell_iter->GetCellData()->SetItem("BranchingCell", 0.0);
        cell_iter->GetCellData()->SetItem("ConsecutiveBranches", 1.0); 
        cell_iter->GetCellData()->SetItem("BirthTimeCell", t);
        cell_iter->GetCellData()->SetItem("TortuosityParent", tortuosity);
        cell_iter->GetCellData()->SetItem("TortuosityDaughter", tortuosity);
        cell_iter->GetCellData()->SetItem("LoopNumber", 0.0);
    }

    UpdateCellData(rCellPopulation);
}


template<unsigned DIM>
void DaughterCellModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    NodeBasedCellPopulation<DIM>* p_node_population = dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);
    if (!p_node_population)
    {
        EXCEPTION("NodeBasedCellPopulation only valid for node");
    }

    // we go through all the cells from the cell population
    for(typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin(); cell_iter != rCellPopulation.End(); ++cell_iter)
    {
        // we need to check if the cell is a daughter or not //

        // we collect the cell data (node index and cell pointer)
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        CellPtr pCell = rCellPopulation.GetCellUsingLocationIndex(node_index);

        // we collect the node indice in CellData (which corresponds to the parent cell indice)
        unsigned node_index_data = pCell->GetCellData()->GetItem("node_index");
        CellPtr pCellData = rCellPopulation.GetCellUsingLocationIndex(node_index_data);

        // we check if it is the same than the one of the cell : if yes, it is the parent cell, if no, it is the daughter cell 
        if(node_index != node_index_data){
            CellPtr pDaughterCell = pCell;
            CellPtr pParentCell = pCellData;
            
            // we collect the daughter type of cell 
            double type_cell = pParentCell->GetCellData()->GetItem("daughter_type");

            if(type_cell == 1.0){
                ++mHighestBranch;

                // the new cell is a tip cell
                MAKE_PTR(TipCellMutationState, p_tip_state);
                pDaughterCell->SetMutationState(p_tip_state);
                MAKE_PTR(StemCellProliferativeType, p_tip_type);
                pDaughterCell->SetCellProliferativeType(p_tip_type);
                pDaughterCell->GetCellData()->SetItem("BranchNumber", mHighestBranch);
                pDaughterCell->GetCellData()->SetItem("BranchingPoint", 0);
                pDaughterCell->GetCellData()->SetItem("LoopNumber", 0.0);

                unsigned branching_cell_parent = pParentCell->GetCellData()->GetItem("BranchingCell");
                pDaughterCell->GetCellData()->SetItem("BranchingCell", branching_cell_parent);

                double new_consecutive_branches_number = pParentCell->GetCellData()->GetItem("ConsecutiveBranches") + 1.0;
                pParentCell->GetCellData()->SetItem("ConsecutiveBranches", new_consecutive_branches_number);
                pDaughterCell->GetCellData()->SetItem("ConsecutiveBranches", new_consecutive_branches_number);

                double t = SimulationTime::Instance()->GetTime(); 
                pDaughterCell->GetCellData()->SetItem("BirthTimeCell", t);
                pParentCell->GetCellData()->SetItem("BirthTimeCell", t);

                ++mHighestBranch;
                pParentCell->GetCellData()->SetItem("BranchNumber", mHighestBranch);

            } else {
                // the new cell is a vessel segment 
                MAKE_PTR(VesselCellMutationState, p_vessel_state);
                pDaughterCell->SetMutationState(p_vessel_state);
                MAKE_PTR(DifferentiatedCellProliferativeType, p_vessel_type);
                pDaughterCell->SetCellProliferativeType(p_vessel_type);
                pDaughterCell->InitialiseCellCycleModel();
                pDaughterCell->GetCellData()->SetItem("BranchingPoint", 0);
            }
            pDaughterCell->GetCellData()->SetItem("node_index", node_index);
            pDaughterCell->GetCellData()->SetItem("OriginalParent", node_index_data);

            double tortuosity = pDaughterCell->GetCellData()->GetItem("TortuosityDaughter");
            pDaughterCell->GetCellData()->SetItem("TortuosityParent", tortuosity);
        }
    }

    // Anastomosis section
    for(typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin(); cell_iter != rCellPopulation.End(); ++cell_iter)
    {
        // we check for possible anastomosis points 
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        CellPtr pCell = rCellPopulation.GetCellUsingLocationIndex(node_index);
        
        if (pCell->GetMutationState()->IsType<TipCellMutationState>()){
            c_vector<double, DIM> cell_position = rCellPopulation.GetLocationOfCellCentre(pCell);
            // anastomosis only possible if cells far enough from the origin and from its branching cell 
            // anastomosis also only possible after a few time step in this position (enough stress applied to the cell)
            CalculateAnastomosisVector(rCellPopulation, p_node_population, pCell);
        }
    }
}

template<unsigned DIM>
void DaughterCellModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class DaughterCellModifier<1>;
template class DaughterCellModifier<2>;
template class DaughterCellModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DaughterCellModifier)

