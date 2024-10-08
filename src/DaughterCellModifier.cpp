#include "DaughterCellModifier.hpp"
#include "MeshBasedCellPopulation.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "CellData.hpp"

#include "TipCellMutationState.hpp"
#include "VesselCellMutationState.hpp"

#include "RandomNumberGenerator.hpp"
#include <algorithm>

#include "Debug.hpp" 
#include "SimulationTime.hpp"

template<unsigned DIM>
DaughterCellModifier<DIM>::DaughterCellModifier(double OldNumNodes)
    : AbstractCellBasedSimulationModifier<DIM>(), mHighestBranch(0)
{
    mOldNumNodes = OldNumNodes;
}

template<unsigned DIM>
DaughterCellModifier<DIM>::~DaughterCellModifier()
{
}

template<unsigned DIM>
double DaughterCellModifier<DIM>::GetOldNumberofNodes(){
    return mOldNumNodes;
}

template<unsigned DIM>
void DaughterCellModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void DaughterCellModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */

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
        cell_iter->GetCellData()->SetItem("daughter_type", 0.0);
        cell_iter->GetCellData()->SetItem("BranchNumber", 0.0);
        cell_iter->GetCellData()->SetItem("BranchingPoint", 0.0);
        cell_iter->GetCellData()->SetItem("DivisionNumber", 0.0);
        cell_iter->GetCellData()->SetItem("BranchingLeader", 0.0);
        cell_iter->GetCellData()->SetItem("ConsecutiveBranches", 1.0);
        cell_iter->GetCellData()->SetItem("BirthTimeCell", t);
        cell_iter->GetCellData()->SetItem("TortuosityParent", tortuosity);
        cell_iter->GetCellData()->SetItem("TortuosityDaughter", tortuosity);
    }

    UpdateCellData(rCellPopulation);
}


template<unsigned DIM>
void DaughterCellModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    //TRACE("Begin Daughter Cell Modifier");

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

    //TRACE("End Daughter Cell Modifier");
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

