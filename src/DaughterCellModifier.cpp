#include "DaughterCellModifier.hpp"
#include "MeshBasedCellPopulation.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "CellData.hpp"

#include "RandomNumberGenerator.hpp"
#include <algorithm>

#include "Debug.hpp" 

template<unsigned DIM>
DaughterCellModifier<DIM>::DaughterCellModifier(double OldNumNodes)
    : AbstractCellBasedSimulationModifier<DIM>()
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
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
double DaughterCellModifier<DIM>::GetAngleFromVectors(c_vector<double,DIM> u, c_vector<double,DIM> v)
{
    // initialisation  
    double scalar_product_uv;
    double angle_uv; 

    // we start by calculating the scalar product between u and v 
    if(DIM == 2){
        scalar_product_uv = u[0]*v[0] + u[1]*v[1];
    } else if(DIM == 3){
        scalar_product_uv = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
    } else{
        scalar_product_uv = u[0]*v[0];
    }

    if(norm_2(u)!= 0 && norm_2(v) != 0){
        double value_uv = scalar_product_uv/(norm_2(u)*norm_2(v));
        // we have to compare the previous value in case it is not between -1 and 1 (due to computational error)
        double hack = std::min(1.0, std::max(-1.0, value_uv));
        angle_uv = std::acos(hack);
        PRINT_4_VARIABLES(angle_uv, scalar_product_uv, norm_2(u), norm_2(v));
    } else {
        angle_uv = 0;
    }

    return angle_uv;
}

template<unsigned DIM>
double DaughterCellModifier<DIM>::ClosestAngleVesselSegment(AbstractCellPopulation<DIM, DIM>& rCellPopulation, CellPtr cell_ptr, std::set<unsigned> neighbouring_node_indices)
{
    c_vector<double,DIM> xj = rCellPopulation.GetLocationOfCellCentre(cell_ptr);

    double alphangular = M_PI/2; // we initialise such that the angle is different from ±PI and 0 : case of a tip cell
    c_vector<double, DIM> u;
    c_vector<double, DIM> v;

    // test : to remove later 
    if (neighbouring_node_indices.size() < 1) {
        TRACE("0 neighbour found");
        PRINT_VARIABLE(neighbouring_node_indices.size());
        EXCEPTION("DAUGHTER CELL MODIFIER");
    }

    if (neighbouring_node_indices.size() > 1) {
        CellPtr p_neighbour_cell_i = rCellPopulation.GetCellUsingLocationIndex(*neighbouring_node_indices.begin());
        CellPtr p_neighbour_cell_k = rCellPopulation.GetCellUsingLocationIndex(*neighbouring_node_indices.end());
        c_vector<double,DIM> xi = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_i);
        c_vector<double,DIM> xk = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_k);
    
        u = rCellPopulation.rGetMesh().GetVectorFromAtoB(xi, xj);
        v = rCellPopulation.rGetMesh().GetVectorFromAtoB(xk, xj);

        alphangular = GetAngleFromVectors(u, v);

        if(neighbouring_node_indices.size() > 2){
            // we look for the two closest neighbours and consider that the angle they do together is alphangular
            c_vector<double, DIM> xmin1 = xi;
            c_vector<double, DIM> xmin2 = xk;
            double length1 = norm_2(rCellPopulation.rGetMesh().GetVectorFromAtoB(xmin1, xj));
            double length2 = norm_2(rCellPopulation.rGetMesh().GetVectorFromAtoB(xmin2, xj));
            double lengthmin1 = std::min(length1, length2);
            double lengthmin2 = std::max(length1, length2);
            
            for(std::set<unsigned>::iterator k = neighbouring_node_indices.begin();
            k != neighbouring_node_indices.end();
            ++k){
                CellPtr p_neighbour_cell_k = rCellPopulation.GetCellUsingLocationIndex(*k);
                c_vector<double, DIM> xk = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_k);
                c_vector<double, DIM> xkj = rCellPopulation.rGetMesh().GetVectorFromAtoB(xk, xj);
                double length_neighbour = norm_2(xkj);
                if(length_neighbour < lengthmin1){
                    lengthmin2 = lengthmin1;
                    xmin2 = xmin1;
                    lengthmin1 = length_neighbour;
                    xmin1 = xk;
                } else if (length_neighbour > lengthmin1 && length_neighbour < lengthmin2){
                    lengthmin2 = length_neighbour;
                    xmin2 = xk;
                }
            }

            // once we have the two closest neighbours then we can calculate the angle alphangular 
            c_vector<double, DIM> u = rCellPopulation.rGetMesh().GetVectorFromAtoB(xmin1, xj);;
            c_vector<double, DIM> v = rCellPopulation.rGetMesh().GetVectorFromAtoB(xmin2, xj);;
            alphangular = GetAngleFromVectors(u, v);
        }
    }

    return alphangular;
}

template<unsigned DIM>
double DaughterCellModifier<DIM>::OptimalAngleVesselSegment(AbstractCellPopulation<DIM, DIM>& rCellPopulation, CellPtr cell_ptr, std::set<unsigned> neighbouring_node_indices)
{
    c_vector<double,DIM> xj = rCellPopulation.GetLocationOfCellCentre(cell_ptr);

    double alphangular = M_PI; // we initialise such as the force is equal to 0, geometrically, it means that the elements are aligned
    c_vector<double, DIM> u;
    c_vector<double, DIM> v;

    CellPtr p_neighbour_cell_i = rCellPopulation.GetCellUsingLocationIndex(*neighbouring_node_indices.begin());
    CellPtr p_neighbour_cell_k = rCellPopulation.GetCellUsingLocationIndex(*neighbouring_node_indices.end());
    c_vector<double,DIM> xi = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_i);
    c_vector<double,DIM> xk = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_k);
    
    u = rCellPopulation.rGetMesh().GetVectorFromAtoB(xi, xj);
    v = rCellPopulation.rGetMesh().GetVectorFromAtoB(xk, xj);

    double alphangularmin = GetAngleFromVectors(u, v);

    if(neighbouring_node_indices.size() > 2){
        // if there are more than two neighbours, then we need to sort out which one are making the smallest angle 
        c_vector<double, DIM> ximin = xi;
        c_vector<double, DIM> xkmin = xk;
        double alpha;
        
        for(std::set<unsigned>::iterator i = neighbouring_node_indices.begin();
        i != neighbouring_node_indices.end();
        ++i){
            for(std::set<unsigned>::iterator k = neighbouring_node_indices.begin();
            k != neighbouring_node_indices.end();
            ++k){
                if(i != k){
                    CellPtr p_neighbour_cell_i = rCellPopulation.GetCellUsingLocationIndex(*i);
                    CellPtr p_neighbour_cell_k = rCellPopulation.GetCellUsingLocationIndex(*k);
                    c_vector<double, DIM> xi = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_i);
                    c_vector<double, DIM> xk = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_k);
                    c_vector<double, DIM> xij = rCellPopulation.rGetMesh().GetVectorFromAtoB(xi, xj);
                    c_vector<double, DIM> xkj = rCellPopulation.rGetMesh().GetVectorFromAtoB(xk, xj);
                    alpha = GetAngleFromVectors(xij, xkj);

                    if(std::abs(alpha) < std::abs(alphangularmin)){
                        alphangularmin = alpha;
                        ximin = xi;
                        xkmin = xk;
                    }
                }
            }
        } 
        u = rCellPopulation.rGetMesh().GetVectorFromAtoB(ximin, xj);
        v = rCellPopulation.rGetMesh().GetVectorFromAtoB(xkmin, xj);
        alphangular = GetAngleFromVectors(u, v);
    } else if(neighbouring_node_indices.size() == 2){
        alphangular = alphangularmin;
        TRACE("test");
    }

    return alphangular;
}

// function that returns the approximation of the length of the vessel element k
template<unsigned DIM>
double DaughterCellModifier<DIM>::LengthVesselSegment(AbstractCellPopulation<DIM, DIM>& rCellPopulation, CellPtr cell_ptr, std::set<unsigned> neighbouring_node_indices){
    double length_vessel_segment = 0; 

    // select the element 
    c_vector<double,DIM> xj = rCellPopulation.GetLocationOfCellCentre(cell_ptr);

    // consider the neighbours of the element 
    if(neighbouring_node_indices.size() > 2){ 
        for(std::set<unsigned>::iterator k = neighbouring_node_indices.begin();
        k != neighbouring_node_indices.end();
        ++k){
            CellPtr p_neighbour_cell_k = rCellPopulation.GetCellUsingLocationIndex(*k);
            c_vector<double, DIM> xk = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_k);

            c_vector<double, DIM> u = rCellPopulation.rGetMesh().GetVectorFromAtoB(xk, xj);
            length_vessel_segment = length_vessel_segment + norm_2(u);
        }
        length_vessel_segment = length_vessel_segment/neighbouring_node_indices.size();
    } else if(neighbouring_node_indices.size() == 2){
        CellPtr p_neighbour_cell_i = rCellPopulation.GetCellUsingLocationIndex(*neighbouring_node_indices.begin());
        CellPtr p_neighbour_cell_k = rCellPopulation.GetCellUsingLocationIndex(*neighbouring_node_indices.end());
        c_vector<double,DIM> xi = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_i);
        c_vector<double,DIM> xk = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_k);

        c_vector<double, DIM> u = rCellPopulation.rGetMesh().GetVectorFromAtoB(xi, xj);
        c_vector<double, DIM> v = rCellPopulation.rGetMesh().GetVectorFromAtoB(xk, xj);
        length_vessel_segment = (norm_2(u)+norm_2(v))/2; 
    }

    // depending on the number of neighbours, we calculate the length of the vessel element
    return length_vessel_segment;
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
            int type_cell = pParentCell->GetCellData()->GetItem("daughter_type");

            if(type_cell == 1.0){
                // the new cell is a tip cell 
                boost::shared_ptr<AbstractCellProperty> p_tipcell_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
                pDaughterCell->SetCellProliferativeType(p_tipcell_type);
                pDaughterCell->GetCellData()->SetItem("BranchLeader", node_index);
                pDaughterCell->GetCellData()->SetItem("BranchingPoint", 0);
                //pParentCell->GetCellData()->SetItem("BranchLeader", -10.0);
            } else {
                // the new cell is a vessel segment 
                boost::shared_ptr<AbstractCellProperty> p_vesselsegment_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());
                pDaughterCell->SetCellProliferativeType(p_vesselsegment_type);
                pDaughterCell->GetCellData()->SetItem("BranchingPoint", 0);
                if(pParentCell->GetCellData()->GetItem("BranchingPoint") == 1.0){
                    unsigned branch_leader_daughter = pParentCell->GetCellData()->GetItem("BranchNeighbourLeader");
                    pDaughterCell->GetCellData()->SetItem("BranchLeader", branch_leader_daughter);
                }
            }
            pDaughterCell->GetCellData()->SetItem("node_index", node_index);
            pDaughterCell->GetCellData()->SetItem("OriginalParent", node_index_data);
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

