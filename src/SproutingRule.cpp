#include "SproutingRule.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "RandomNumberGenerator.hpp"
#include <algorithm>
#include "Debug.hpp"
#include <random>
#include <cstdlib>
#include <iterator>


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SproutingRule<ELEMENT_DIM, SPACE_DIM>::SproutingRule(double Psprout, double directioncoeff, double perpendiculardirectioncoeff)
    : AbstractCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>()
{
    assert(Psprout > 0 && Psprout < 1);
    mPsprout = Psprout;
    mDirectionCoeff = directioncoeff;
    mPerpendicularDirectionCoeff = perpendiculardirectioncoeff;
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
double SproutingRule<ELEMENT_DIM, SPACE_DIM>::GetDirectionCoefficient(){
    return mDirectionCoeff;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double SproutingRule<ELEMENT_DIM, SPACE_DIM>::GetPerpendicularDirectionCoefficient(){
    return mPerpendicularDirectionCoeff;
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

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::tuple<CellPtr, unsigned, c_vector<double,SPACE_DIM>> SproutingRule<ELEMENT_DIM, SPACE_DIM>::RandomNeighbour(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pCell, std::set<unsigned> neighbouring_node_indices)
{
    // we collect the cell data (node index)
    unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(pCell);
    double branching_point = pCell->GetCellData()->GetItem("BranchingPoint");

    // we start by creating a new set with the cell's neighbour from the same branch 
    std::set<unsigned> neighbouring_node_indices_same_branch;
    for(std::set<unsigned>::iterator k = neighbouring_node_indices.begin(); k != neighbouring_node_indices.end(); ++k){
        // we collect the cell pointer 
        CellPtr pNeighbourCell = rCellPopulation.GetCellUsingLocationIndex(*k);
        unsigned original_parent_neighbour = static_cast<unsigned>(pNeighbourCell->GetCellData()->GetItem("OriginalParent"));
        if (pNeighbourCell->GetCellData()->GetItem("BranchLeader") == pCell->GetCellData()->GetItem("BranchLeader")){
            // we add the cell to the new neighbour set 
            neighbouring_node_indices_same_branch.insert(*k);
        } else if (branching_point == 1.0 && original_parent_neighbour == node_index){
            neighbouring_node_indices_same_branch.insert(*k);
        }
    }

    // we generate a random number to obtain a random neighbour indice 
    int random_index = std::rand() % neighbouring_node_indices_same_branch.size();

    // we collect the neighbour indice corresponding to the random index generated 
    unsigned random_indice = GetNthNeighbourIndice(neighbouring_node_indices_same_branch, random_index);

    // we collect the data from the random cell (cell pointer and coordinates)
    CellPtr pRandomNeighbour = rCellPopulation.GetCellUsingLocationIndex(random_indice);
    c_vector<double, SPACE_DIM> random_neighbour = rCellPopulation.GetLocationOfCellCentre(pRandomNeighbour);

    // depending on the number of neighbours, we calculate the length of the vessel element
    return std::make_tuple(pRandomNeighbour, random_indice, random_neighbour);
}

// function that returns the closest neighbour of the vessel element k
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> SproutingRule<ELEMENT_DIM, SPACE_DIM>::ClosestNeighbour(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pCell, std::set<unsigned> neighbouring_node_indices){
    // we collect the cell's neighbours data (pointers, coordinates)
    CellPtr p_closest_neighbour = rCellPopulation.GetCellUsingLocationIndex(*neighbouring_node_indices.begin());
    c_vector<double, SPACE_DIM> closest_neighbour = rCellPopulation.GetLocationOfCellCentre(p_closest_neighbour);
    double length_closest_neighbour = norm_2(closest_neighbour); 

    // select the element 
    c_vector<double,SPACE_DIM> xj = rCellPopulation.GetLocationOfCellCentre(pCell);

    // consider the neighbours of the element 
    for(std::set<unsigned>::iterator k = neighbouring_node_indices.begin();
    k != neighbouring_node_indices.end();
    ++k){
        CellPtr p_neighbour_cell_k = rCellPopulation.GetCellUsingLocationIndex(*k);
        c_vector<double, SPACE_DIM> xk = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_k);

        c_vector<double, SPACE_DIM> u = rCellPopulation.rGetMesh().GetVectorFromAtoB(xk, xj);
        double length_u = norm_2(u);
        
        if(length_u < length_closest_neighbour){
            length_closest_neighbour = length_u;
            closest_neighbour = xk;
        }
    }

    // depending on the number of neighbours, we calculate the length of the vessel element
    return closest_neighbour;
}

// function that returns the closest neighbour of the vessel element k
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> SproutingRule<ELEMENT_DIM, SPACE_DIM>::TwoClosestNeighbours(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pCell, std::set<unsigned> neighbouring_node_indices){
    // test : to be removed later 
    TRACE("Begin TwoClosestNeighbours");
    PRINT_VARIABLE(neighbouring_node_indices.size());

    unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(pCell);
    PRINT_VARIABLE(node_index);

    // initialisation 
    c_vector<double, SPACE_DIM> closest_neighbours_vector;

    c_vector<double, SPACE_DIM> xj = rCellPopulation.GetLocationOfCellCentre(pCell);

    // we need to check the number of neighbours 
    if(neighbouring_node_indices.size() == 2){
        // we collect the coordinates of the 2 first node neighbours 
        CellPtr p_neighbour_cell_1 = rCellPopulation.GetCellUsingLocationIndex(*neighbouring_node_indices.begin());
        unsigned second_neighbour_node_indice = GetNthNeighbourIndice(neighbouring_node_indices, 1);
        CellPtr p_neighbour_cell_2 = rCellPopulation.GetCellUsingLocationIndex(second_neighbour_node_indice);
        c_vector<double, SPACE_DIM> x1 = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_1);
        c_vector<double, SPACE_DIM> x2 = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_2);

        closest_neighbours_vector = rCellPopulation.rGetMesh().GetVectorFromAtoB(x1, x2);

        // test : to be removed later 
        PRINT_VECTOR(xj);
        PRINT_VECTOR(x1);
        PRINT_VECTOR(x2);
        PRINT_VECTOR(closest_neighbours_vector);
    
    } else if (neighbouring_node_indices.size() > 2) {
        // we collect the coordinates of the 2 first node neighbours 
        CellPtr p_neighbour_cell_i = rCellPopulation.GetCellUsingLocationIndex(*neighbouring_node_indices.begin());
        unsigned second_neighbour_node_indice = GetNthNeighbourIndice(neighbouring_node_indices, 1);
        CellPtr p_neighbour_cell_k = rCellPopulation.GetCellUsingLocationIndex(second_neighbour_node_indice);
        c_vector<double, SPACE_DIM> xi = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_i);
        c_vector<double, SPACE_DIM> xk = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_k);

        // we look for the two closest neighbours and consider that the angle they do together is alphangular
        c_vector<double, SPACE_DIM> xmin1 = xi;
        c_vector<double, SPACE_DIM> xmin2 = xk;
        double length1 = norm_2(rCellPopulation.rGetMesh().GetVectorFromAtoB(xmin1, xj));
        double length2 = norm_2(rCellPopulation.rGetMesh().GetVectorFromAtoB(xmin2, xj));
        double lengthmin1 = std::min(length1, length2);
        double lengthmin2 = std::max(length1, length2);
        
        for(std::set<unsigned>::iterator k = neighbouring_node_indices.begin(); k != neighbouring_node_indices.end(); ++k){
            CellPtr p_neighbour_cell_k = rCellPopulation.GetCellUsingLocationIndex(*k);
            c_vector<double, SPACE_DIM> xk = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_k);
            c_vector<double, SPACE_DIM> xkj = rCellPopulation.rGetMesh().GetVectorFromAtoB(xk, xj);
            double length_neighbour = norm_2(xkj);
            if(length_neighbour > 0 && length_neighbour < lengthmin1){
                lengthmin2 = lengthmin1;
                xmin2 = xmin1;
                lengthmin1 = length_neighbour;
                xmin1 = xk;
            } else if (length_neighbour > 0 && length_neighbour > lengthmin1 && length_neighbour < lengthmin2){
                lengthmin2 = length_neighbour;
                xmin2 = xk;
            }
        }

        // once we have the two closest neighbours then we can calculate the vector between them 
        closest_neighbours_vector = rCellPopulation.rGetMesh().GetVectorFromAtoB(xmin1, xmin2);

        // test : to be removed later 
        PRINT_VECTOR(xj);
        PRINT_VECTOR(xi);
        PRINT_VECTOR(xk);
        PRINT_VECTOR(xmin1);
        PRINT_VECTOR(xmin2);
        PRINT_VECTOR(closest_neighbours_vector);

    } else if (neighbouring_node_indices.size() < 2){
        closest_neighbours_vector = ClosestNeighbour(rCellPopulation, pCell, neighbouring_node_indices);
    }

    // test : to be removed later 
    TRACE("End TwoClosestNeighbours");

    // depending on the number of neighbours, we calculate the length of the vessel element
    // we return the vector between the two closest neighbours  
    return closest_neighbours_vector;
}

// Overridden CalculateCellDivisionVector() method, return the daughter cell positions
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > SproutingRule<ELEMENT_DIM, SPACE_DIM>::CalculateCellDivisionVector(CellPtr pParentCell, AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    // test : to be removed 
    TRACE("Begin Division Rule");

    // initialisation 
    c_vector<double, SPACE_DIM> parent_position;
    c_vector<double, SPACE_DIM> daughter_position;

    // consider a parent cell 
    //unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(pParentCell);
    parent_position = rCellPopulation.GetLocationOfCellCentre(pParentCell);
    daughter_position = zero_vector<double>(SPACE_DIM);

    // we generate a probability for the cell to be a tip cell and create a new sprout 
    // the location of the daughter cell will depend on the type of the daughter cell 
    // if the daughter cell is a vessel segment, we incorporate it directly inside the vessel 
    // if the daughter cell is a tip cell, we add it perpendicularly to the vector between the parent cell and its closest neighbour 

    // we generate a random number that will set the type of daughter cell 
    double random_proba = RandomNumberGenerator::Instance()->ranf();

    PRINT_VARIABLE(random_proba);

    // we check the neighbours of the parent cell 
    NodeBasedCellPopulation<SPACE_DIM>* p_node_population = dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation);
    if (!p_node_population)
    {
        EXCEPTION("NodeBasedCellPopulation only valid for node");
    }
    std::set<unsigned> neighbouring_node_indices = p_node_population->GetNeighbouringLocationIndices(pParentCell);

    // test : to be removed 
    PRINT_VARIABLE(neighbouring_node_indices.size());

    // we look for a random neighbour to have a landmark 
    //c_vector<double, SPACE_DIM> closest_neighbour = ClosestNeighbour(rCellPopulation, pParentCell, neighbouring_node_indices);
    CellPtr pRandomNeighbour;
    unsigned random_indice;
    c_vector<double, SPACE_DIM> random_neighbour;
    std::tie(pRandomNeighbour, random_indice, random_neighbour) = RandomNeighbour(rCellPopulation, pParentCell, neighbouring_node_indices);
    unsigned branch_leader_random_neighbour_indice = static_cast<unsigned>(pRandomNeighbour->GetCellData()->GetItem("BranchLeader"));
    CellPtr p_branch_leader_random_neighbour = rCellPopulation.GetCellUsingLocationIndex(branch_leader_random_neighbour_indice);
    c_vector<double, SPACE_DIM> branch_leader_random_neighbour_position = rCellPopulation.GetLocationOfCellCentre(p_branch_leader_random_neighbour); 

    // if the daughter cell is a vessel segment, its position will be inside the vessel, hence, along the vector between the parent cell and its closest neighbour 
    //c_vector<double, SPACE_DIM> daughter_direction = rCellPopulation.rGetMesh().GetVectorFromAtoB(parent_position, branch_leader_random_neighbour_position)/norm_2(rCellPopulation.rGetMesh().GetVectorFromAtoB(parent_position, branch_leader_random_neighbour_position));
    c_vector<double, SPACE_DIM> daughter_direction = rCellPopulation.rGetMesh().GetVectorFromAtoB(parent_position, random_neighbour)/norm_2(rCellPopulation.rGetMesh().GetVectorFromAtoB(parent_position, random_neighbour));

    // if the daughter cell is a tip cell, its position will be perpendicular to the vector between the parent cell and its closest neighbour 
    c_vector<double, SPACE_DIM> perpendicular_daughter_direction = zero_vector<double>(SPACE_DIM);

    // test : to be removed 
    PRINT_VECTOR(daughter_direction);
    
    if(SPACE_DIM == 2){
        // counterclockwise rotation 
        perpendicular_daughter_direction(0) = -daughter_direction(1);
        perpendicular_daughter_direction(1) = daughter_direction(0);
    } else if(SPACE_DIM == 3){
        // counterclockwise rotation around x-axis 
        perpendicular_daughter_direction(0) = daughter_direction(0);
        perpendicular_daughter_direction(1) = -daughter_direction(2);
        perpendicular_daughter_direction(2) = daughter_direction(1);
    } else {
        perpendicular_daughter_direction(0) = -daughter_direction(0);
    }

    if(random_proba < mPsprout){
        // tip cell : new sprout 
        // addition of the daughter cell to the map
        daughter_position = parent_position + mPerpendicularDirectionCoeff*perpendicular_daughter_direction;
        pParentCell->GetCellData()->SetItem("daughter_type", 1.0);
        pParentCell->GetCellData()->SetItem("BranchingPoint", 1.0);
    } else {
        // vessel segment : incorporation in the vessel 
        // addition of the daughter cell to the map
        daughter_position = parent_position + mDirectionCoeff*daughter_direction;
        pParentCell->GetCellData()->SetItem("daughter_type", 0.0);
        if(pParentCell->GetCellData()->GetItem("BranchingPoint") == 1.0){
            pParentCell->GetCellData()->SetItem("BranchNeighbourLeader", random_indice);
        }
    }

    double division_number = pParentCell->GetCellData()->GetItem("DivisionNumber");
    double new_division_number = division_number + 1.0;
    pParentCell->GetCellData()->SetItem("DivisionNumber", new_division_number);

    // test : to be removed 
    PRINT_VECTOR(parent_position);
    PRINT_VECTOR(daughter_position);

    // we return the coordinates of the parent cell and daughter cell 
    std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > positions(parent_position, daughter_position);

    // test : to be removed 
    TRACE("End Division Rule");

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