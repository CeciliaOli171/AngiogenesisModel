#include "SproutingRule.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "VesselCellMutationState.hpp"
#include "TipCellMutationState.hpp"
#include "BranchingCellMutationState.hpp"
#include "UniformCellCycleModel.hpp"

#include "RandomNumberGenerator.hpp"
#include <algorithm>
#include "Debug.hpp"
#include <random>
#include <cstdlib>
#include <iterator>


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SproutingRule<ELEMENT_DIM, SPACE_DIM>::SproutingRule(double Psprout)
    : AbstractCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>()
{
    assert(Psprout > 0 && Psprout < 1);
    mPsprout = Psprout;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SproutingRule<ELEMENT_DIM, SPACE_DIM>::~SproutingRule()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double SproutingRule<ELEMENT_DIM, SPACE_DIM>::GetSproutingProbability(){
    return mPsprout;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> SproutingRule<ELEMENT_DIM, SPACE_DIM>::GetVesselSegmentNeighbours(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation,NodeBasedCellPopulation<SPACE_DIM>* p_node_population, CellPtr pParentCell){  
    std::set<unsigned> neighbours_set = p_node_population->GetNodesWithinNeighbourhoodRadius(pParentCell->GetCellId(),1.5);
    std::set<unsigned> neighbours_set_without_tip;

    for(std::set<unsigned>::iterator k = neighbours_set.begin(); k != neighbours_set.end(); ++k){
        // we collect the cell pointer 
        CellPtr pNeighbourCell = rCellPopulation.GetCellUsingLocationIndex(*k);
        if (pNeighbourCell->GetMutationState()->IsType<TipCellMutationState>() == false){
            // we add the cell to the new neighbour set 
            neighbours_set_without_tip.insert(*k);
        } 
    }
    return neighbours_set_without_tip;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned SproutingRule<ELEMENT_DIM, SPACE_DIM>::GetNthNeighbourIndice(std::set<unsigned> neighbouring_node_indices, int n){
    // initialisation
    int neighbour_indice = *neighbouring_node_indices.begin();

    int iter = 0;
    for(std::set<unsigned>::iterator k = neighbouring_node_indices.begin();
    k != neighbouring_node_indices.end();
    ++k){
        if(iter == n){
            neighbour_indice = *k;
        }
        ++iter;
    }

    return neighbour_indice;
}

// function that returns the closest neighbour of the vessel element k
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std:: pair <c_vector<double, SPACE_DIM>, unsigned> SproutingRule<ELEMENT_DIM, SPACE_DIM>::ClosestNeighbour(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pCell, std::set<unsigned> neighbouring_node_indices){
    // we collect the cell's neighbours data (pointers, coordinates)
    CellPtr p_closest_neighbour = rCellPopulation.GetCellUsingLocationIndex(*neighbouring_node_indices.begin());
    c_vector<double, SPACE_DIM> closest_neighbour = rCellPopulation.GetLocationOfCellCentre(p_closest_neighbour);
    double length_closest_neighbour = norm_2(closest_neighbour); 

    // select the element 
    c_vector<double,SPACE_DIM> xj = rCellPopulation.GetLocationOfCellCentre(pCell);

    unsigned indice_closest_neighbour = *neighbouring_node_indices.begin();

    // consider the neighbours of the element 
    for(std::set<unsigned>::iterator k = neighbouring_node_indices.begin();
    k != neighbouring_node_indices.end();
    ++k){
        CellPtr p_neighbour_cell_k = rCellPopulation.GetCellUsingLocationIndex(*k);
        c_vector<double, SPACE_DIM> xk = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_k);

        c_vector<double, SPACE_DIM> u = rCellPopulation.rGetMesh().GetVectorFromAtoB(xk, xj);
        double length_u = norm_2(u);

        // check if it is a branching point or in the same branch and remove tip cells 
        if(p_neighbour_cell_k->GetCellData()->GetItem("BranchNumber") == pCell->GetCellData()->GetItem("BranchNumber")){
            if(length_u < length_closest_neighbour){
            length_closest_neighbour = length_u;
            closest_neighbour = xk;
            indice_closest_neighbour = *k;
            }
        } else if(p_neighbour_cell_k->GetMutationState()->IsType<BranchingCellMutationState>() && pCell->GetCellData()->GetItem("BranchingLeader") == *k){
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


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > SproutingRule<ELEMENT_DIM, SPACE_DIM>::CalculateCellDivisionVesselCellVector(CellPtr pParentCell, AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    // initialisation 
    c_vector<double, SPACE_DIM> parent_position;
    c_vector<double, SPACE_DIM> daughter_position;

    // consider a parent cell 
    parent_position = rCellPopulation.GetLocationOfCellCentre(pParentCell);
    daughter_position = zero_vector<double>(SPACE_DIM);

    // we return the coordinates of the parent cell and daughter cell 
    std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > positions(parent_position, daughter_position);

    TRACE("Test NO division of vessel segment");

    return positions;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> SproutingRule<ELEMENT_DIM, SPACE_DIM>::PerpendicularDaughterDirection(c_vector<double, SPACE_DIM> daughter_direction, c_vector<double, SPACE_DIM> parent_position){
    // initialisation of the perpendicular vector 
    c_vector<double, SPACE_DIM> perpendicular_daughter_direction;

    if(SPACE_DIM == 2){
        // counterclockwise rotation 
        perpendicular_daughter_direction(0) = -daughter_direction(1);
        perpendicular_daughter_direction(1) = daughter_direction(0);
    } else if(SPACE_DIM == 3){
        // random choice of a, b, c in the interval [-1, 1]
        double a = -1 + 2*RandomNumberGenerator::Instance()->ranf();
        double b = -1 + 2*RandomNumberGenerator::Instance()->ranf();
        double c = -1 + 2*RandomNumberGenerator::Instance()->ranf();

        // orthogonal vector 
        perpendicular_daughter_direction(0) = a*(-daughter_direction(1)) + c*(-daughter_direction(2));
        perpendicular_daughter_direction(1) = a*daughter_direction(0) + b*(-daughter_direction(0));
        perpendicular_daughter_direction(2) = b*daughter_direction(1) + c*daughter_direction(0);
    } else {
        perpendicular_daughter_direction(0) = -daughter_direction(0);
    }

    return perpendicular_daughter_direction;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > SproutingRule<ELEMENT_DIM, SPACE_DIM>::CalculateGrowthVector( AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, NodeBasedCellPopulation<SPACE_DIM>* p_node_population, CellPtr pParentCell){
    // initialisation 
    c_vector<double, SPACE_DIM> parent_position;
    c_vector<double, SPACE_DIM> daughter_position;
    double mDirectionCoeff = rCellPopulation.GetMeinekeDivisionSeparation();
    double mPerpendicularDirectionCoeff = rCellPopulation.GetMeinekeDivisionSeparation();

    // we generate a random number that will set the type of daughter cell 
    double random_proba = RandomNumberGenerator::Instance()->ranf();

    // consider a parent cell 
    parent_position = rCellPopulation.GetLocationOfCellCentre(pParentCell);
    daughter_position = zero_vector<double>(SPACE_DIM);

    // we select the neighbours of the cell that are just vessel segments
    std::set<unsigned> neighbours_set_without_tip = GetVesselSegmentNeighbours(rCellPopulation, p_node_population, pParentCell);

    std::vector<unsigned> neighbouring_node_indices = std::vector<unsigned>(neighbours_set_without_tip.begin(), neighbours_set_without_tip.end());
    
    // we look for a random neighbour to have a landmark 
    if(neighbouring_node_indices.empty()) {
        EXCEPTION("No neighbours. Boo.");
    }

    // we obtain the closest neighbour to the leading tip cell of the branch with the function ClosestNeighbour with the new set of neighbours obtained just before 
    unsigned closest_indice;
    if(neighbouring_node_indices.size() > 1){
        c_vector<double, SPACE_DIM> closest_neighbour;
        std::tie(closest_neighbour, closest_indice) = ClosestNeighbour(rCellPopulation, pParentCell, neighbours_set_without_tip);
    } else if (neighbouring_node_indices.size() == 1){
        closest_indice = *neighbouring_node_indices.begin(); // because it is supposed to have only one neighbour in the same branch 
    }

    // we collect the closest neighbour data (cell pointer, coordinates, indice)
    CellPtr pClosestNeighbour = rCellPopulation.GetCellUsingLocationIndex(closest_indice);
    c_vector<double, SPACE_DIM> closest_neighbour = rCellPopulation.GetLocationOfCellCentre(pClosestNeighbour);

    // if the daughter cell is a vessel segment, its position will be inside the vessel, hence, along the vector between the parent cell and its closest neighbour 
    c_vector<double, SPACE_DIM> daughter_direction = rCellPopulation.rGetMesh().GetVectorFromAtoB(parent_position, closest_neighbour)/norm_2(rCellPopulation.rGetMesh().GetVectorFromAtoB(parent_position, closest_neighbour));
    
    if(random_proba < mPsprout){
        // tip cell : new sprout

        // if the daughter cell is a tip cell, its position will be perpendicular to the vector between the parent cell and its closest neighbour 
        c_vector<double, SPACE_DIM> perpendicular_daughter_direction = PerpendicularDaughterDirection(daughter_direction, parent_position);

        // addition of the daughter cell to the map
        daughter_position = parent_position - 0.5*mPerpendicularDirectionCoeff*perpendicular_daughter_direction + 0.5*mDirectionCoeff*daughter_direction;
        parent_position = parent_position + 0.5*mPerpendicularDirectionCoeff*perpendicular_daughter_direction + 0.5*mDirectionCoeff*daughter_direction;
        pParentCell->GetCellData()->SetItem("daughter_type", 1.0);

        MAKE_PTR(BranchingCellMutationState, p_branching_state);
        pClosestNeighbour->SetMutationState(p_branching_state);
        pClosestNeighbour->GetCellData()->SetItem("BranchingPoint", 1.0);
        pParentCell->GetCellData()->SetItem("BranchingLeader", closest_indice);
        pParentCell->GetCellData()->SetItem("BranchingCell", closest_indice);

        double tortuosity = pParentCell->GetCellData()->GetItem("TortuosityParent");
        double tortuosity_parent = tortuosity + norm_2(parent_position-closest_neighbour);
        double tortuosity_daughter = tortuosity + norm_2(daughter_position-closest_neighbour);
        pParentCell->GetCellData()->SetItem("TortuosityParent", tortuosity_parent);
        pParentCell->GetCellData()->SetItem("TortuosityDaughter", tortuosity_daughter);

        TRACE("NEW SPROUT INCOMING");

        // PRINT_VECTOR(daughter_position);
        // PRINT_VECTOR(parent_position);
        // PRINT_VARIABLE(norm_2(daughter_position-closest_neighbour));
        // PRINT_VARIABLE(norm_2(parent_position-closest_neighbour));

    } else {
        // vessel segment : incorporation in the vessel 
        // addition of the daughter cell to the map
        daughter_position = parent_position;
        parent_position = parent_position - mDirectionCoeff*daughter_direction;
        pParentCell->GetCellData()->SetItem("daughter_type", 0.0);

        double tortuosity = pParentCell->GetCellData()->GetItem("TortuosityParent");
        double tortuosity_parent = tortuosity + norm_2(parent_position-daughter_position);
        double tortuosity_daughter = tortuosity;
        pParentCell->GetCellData()->SetItem("TortuosityParent", tortuosity_parent);
        pParentCell->GetCellData()->SetItem("TortuosityDaughter", tortuosity_daughter);

    }

    double division_number = pParentCell->GetCellData()->GetItem("DivisionNumber");
    double new_division_number = division_number + 1.0;
    pParentCell->GetCellData()->SetItem("DivisionNumber", new_division_number);

    // we return the coordinates of the parent cell and daughter cell 
    std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > positions(parent_position, daughter_position);

    return positions;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > SproutingRule<ELEMENT_DIM, SPACE_DIM>::CalculateCellDivisionTipCellVector(CellPtr pParentCell, AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    // we generate a probability for the cell to be a tip cell and create a new sprout 
    // the location of the daughter cell will depend on the type of the daughter cell 
    // if the daughter cell is a vessel segment, we incorporate it directly inside the vessel 
    // if the daughter cell is a tip cell, we add it perpendicularly to the vector between the parent cell and its closest neighbour 
    std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > positions;

    // we check the neighbours of the parent cell 
    NodeBasedCellPopulation<SPACE_DIM>* p_node_population = dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation);
    if (!p_node_population)
    {
        EXCEPTION("NodeBasedCellPopulation only valid for node");
    }

    positions = CalculateGrowthVector(rCellPopulation, p_node_population, pParentCell);

    return positions;
}

// Overridden CalculateCellDivisionVector() method, return the daughter cell positions
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > SproutingRule<ELEMENT_DIM, SPACE_DIM>::CalculateCellDivisionVector(CellPtr pParentCell, AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    // TRACE("Begin Sprouting Rule"); 

    // initialisation 
    std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > positions;

    // test : to be removed 
    unsigned parent_indice = rCellPopulation.GetLocationIndexUsingCell(pParentCell);

    if (pParentCell->GetMutationState()->IsType<TipCellMutationState>()){
        positions = CalculateCellDivisionTipCellVector(pParentCell, rCellPopulation);
    } else if (pParentCell->GetMutationState()->IsType<VesselCellMutationState>()){
        positions = CalculateCellDivisionVesselCellVector(pParentCell, rCellPopulation);
    } else {
        // case that should never been reached since the branching point cannot divide (differentiated cells)

        // test : to be removed 
        PRINT_VARIABLE(pParentCell->GetAge());
        UniformCellCycleModel* d_cellcycle = dynamic_cast<UniformCellCycleModel*> (pParentCell->GetCellCycleModel());
        PRINT_VARIABLE(d_cellcycle->GetCellCycleDuration());
        PRINT_VARIABLE(pParentCell->GetMutationState()->GetColour());
        PRINT_VARIABLE(pParentCell->GetCellProliferativeType()->GetColour());
        
        // if this case is reached then we send back an error 
        NEVER_REACHED;
    }

    // TRACE("End Sprouting Rule");

    return positions;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SproutingRule<ELEMENT_DIM, SPACE_DIM>::OutputCellCentreBasedDivisionRuleInfo(out_stream& rParamsFile){
    std::string cell_division_rule_type = "Sprouting Division Rule";

    *rParamsFile << "\t\t\t<" << cell_division_rule_type << ">\n";
}


// Explicit instantiation
template class SproutingRule<1,1>;
template class SproutingRule<1,2>;
template class SproutingRule<1,3>;
template class SproutingRule<2,2>;
template class SproutingRule<2,3>;
template class SproutingRule<3,3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SproutingRule)