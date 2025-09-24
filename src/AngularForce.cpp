#include "AngularForce.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include <tuple>

#include "BranchingSegmentMutationState.hpp"
#include "VesselSegmentMutationState.hpp"
#include "VesselTipMutationState.hpp"

#include "Debug.hpp"

template<unsigned DIM>
AngularForce<DIM>::AngularForce(double omegaa)
    : AbstractForce<DIM>()
{
    assert(omegaa >0);
    mOmegaa = omegaa;
}

template<unsigned DIM>
AngularForce<DIM>::~AngularForce()
{
}

template<unsigned DIM>
double AngularForce<DIM>::GetAngularPersistence()
{
    return mOmegaa;
}

template<unsigned DIM>
unsigned AngularForce<DIM>::GetNthNeighbourIndice(std::set<unsigned> neighbouring_node_indices, int n){
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

template<unsigned DIM>
std::set<unsigned> AngularForce<DIM>::GetSameBranchNeighbours(AbstractCellPopulation<DIM, DIM>& rCellPopulation,NodeBasedCellPopulation<DIM>* p_node_population, CellPtr pCell){  
    std::set<unsigned> neighbours_set = p_node_population->GetNodesWithinNeighbourhoodRadius(pCell->GetCellId(),1.5); // we use here the cut-off length between two cells 
    std::set<unsigned> neighbours_set_same_branch;

    for(std::set<unsigned>::iterator k = neighbours_set.begin(); k != neighbours_set.end(); ++k){
        // we collect the cell pointer 
        CellPtr pNeighbourCell = rCellPopulation.GetCellUsingLocationIndex(*k);
        if (pNeighbourCell->GetCellData()->GetItem("BranchNumber") == pCell->GetCellData()->GetItem("BranchNumber")){
            // we add the cell to the new neighbour set 
            neighbours_set_same_branch.insert(*k);
        } 
    }
    return neighbours_set_same_branch;
}

template<unsigned DIM>
double AngularForce<DIM>::GetAngleFromVectors(c_vector<double,DIM> u, c_vector<double,DIM> v)
{
    // initialisation  
    double scalar_product_uv;
    double angle_uv = 0; 

    // we start by calculating the scalar product between u and v 
    for(unsigned i=0; i<DIM; i++){
        scalar_product_uv = u[i]*v[i];
    }

    // we check than none of the norms are equal to 0, if so, we set the angle to be equal to 0
    assert(norm_2(u)*norm_2(v) != 0);
    double value_uv = scalar_product_uv/(norm_2(u)*norm_2(v));

    // we have to compare the previous value in case it is not between -1 and 1 (due to computational error)
    double test_value_uv = std::min(1.0, std::max(-1.0, value_uv));
    angle_uv = std::acos(test_value_uv);

    return angle_uv;
}

template<unsigned DIM>
std::tuple<double, c_vector<double,DIM>, c_vector<double,DIM>> AngularForce<DIM>::ClosestAngleVesselSegment(AbstractCellPopulation<DIM>& rCellPopulation, CellPtr pCell, std::set<unsigned> neighbouring_node_indices)
{
    // initialisation
    c_vector<double, DIM> u;
    c_vector<double, DIM> v;

    double alphangular = M_PI; // we initialise such that the force is equal to 0 if none of the condition are respected 
    
    // coordinates of the cell that we consider 
    c_vector<double,DIM> xj = rCellPopulation.GetLocationOfCellCentre(pCell);

    // we start by creating a new set with the cell's neighbour from the same branch or loop
    std::set<unsigned> neighbouring_node_indices_same_branch;

    for(std::set<unsigned>::iterator k = neighbouring_node_indices.begin(); k != neighbouring_node_indices.end(); ++k){
        // we collect the cell pointer 
        CellPtr pNeighbourCell = rCellPopulation.GetCellUsingLocationIndex(*k);

        if (pNeighbourCell->GetCellData()->GetItem("BranchNumber") == pCell->GetCellData()->GetItem("BranchNumber")){
            // we add the cell to the new neighbour set 
            neighbouring_node_indices_same_branch.insert(*k);
        } else if (pNeighbourCell->GetMutationState() ->IsType<BranchingSegmentMutationState>() && pCell->GetCellData()->GetItem("BranchingSegment") == *k){
            neighbouring_node_indices_same_branch.insert(*k);
        } else if (pNeighbourCell->GetCellData()->GetItem("LoopNumber") == pCell->GetCellData()->GetItem("LoopNumber") && pNeighbourCell->GetCellData()->GetItem("LoopNumber") != 0) {
            // case of the cells in the same loop but with different branches number: use loop number 
            neighbouring_node_indices_same_branch.insert(*k);
        }
    }

    // we need to check the number of neighbours : 
    // if only two neighbours then we consider all the cells even the tip cells 
    // if there are more than two neighbours, then we want to align to the vessel segments
    if(neighbouring_node_indices_same_branch.size() == 2){
        // we collect the coordinates of the 2 first node neighbours 
        CellPtr p_neighbour_cell_i = rCellPopulation.GetCellUsingLocationIndex(*neighbouring_node_indices_same_branch.begin());
        unsigned second_neighbour_indice = GetNthNeighbourIndice(neighbouring_node_indices_same_branch, 1);
        CellPtr p_neighbour_cell_k = rCellPopulation.GetCellUsingLocationIndex(second_neighbour_indice);
        c_vector<double,DIM> xi = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_i);
        c_vector<double,DIM> xk = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_k);
    
        u = rCellPopulation.rGetMesh().GetVectorFromAtoB(xi, xj);
        v = rCellPopulation.rGetMesh().GetVectorFromAtoB(xk, xj);

        alphangular = GetAngleFromVectors(u, v);
    } else if (neighbouring_node_indices_same_branch.size() > 2) {
        // we collect the coordinates of the 2 first node neighbours 
        CellPtr p_neighbour_cell_i = rCellPopulation.GetCellUsingLocationIndex(*neighbouring_node_indices_same_branch.begin());
        unsigned second_neighbour_indice = GetNthNeighbourIndice(neighbouring_node_indices_same_branch, 1);
        CellPtr p_neighbour_cell_k = rCellPopulation.GetCellUsingLocationIndex(second_neighbour_indice);
        c_vector<double,DIM> xi = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_i);
        c_vector<double,DIM> xk = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_k);
    
        u = rCellPopulation.rGetMesh().GetVectorFromAtoB(xi, xj);
        v = rCellPopulation.rGetMesh().GetVectorFromAtoB(xk, xj);

        alphangular = GetAngleFromVectors(u, v);

        // we look for the two closest neighbours and consider that the angle they do together is alphangular
        c_vector<double, DIM> xmin1 = xi;
        c_vector<double, DIM> xmin2 = xk;
        double length1 = norm_2(rCellPopulation.rGetMesh().GetVectorFromAtoB(xmin1, xj));
        double length2 = norm_2(rCellPopulation.rGetMesh().GetVectorFromAtoB(xmin2, xj));
        double lengthmin1 = std::min(length1, length2);
        double lengthmin2 = std::max(length1, length2);
        
        for(std::set<unsigned>::iterator k = neighbouring_node_indices_same_branch.begin();
        k != neighbouring_node_indices_same_branch.end();
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
        c_vector<double, DIM> u = rCellPopulation.rGetMesh().GetVectorFromAtoB(xmin1, xj);
        c_vector<double, DIM> v = rCellPopulation.rGetMesh().GetVectorFromAtoB(xmin2, xj);
        alphangular = GetAngleFromVectors(u, v);
    }

    return std::make_tuple(alphangular, u, v);
}

template<unsigned DIM>
std::tuple<double, c_vector<double,DIM>, c_vector<double,DIM>> AngularForce<DIM>::OptimalAngleVesselElement(AbstractCellPopulation<DIM>& rCellPopulation, CellPtr pCell, std::set<unsigned> neighbouring_node_indices)
{
    // initialisation 
    c_vector<double, DIM> u;
    c_vector<double, DIM> v;

    double alphangular = M_PI; // we initialise such as the force is equal to 0, geometrically, it means that the elements are aligned

    // coordinates of the cell that we consider 
    c_vector<double,DIM> xj = rCellPopulation.GetLocationOfCellCentre(pCell);

    if(neighbouring_node_indices.size() > 1){
        // we consider the two first neighbours of the cell
        CellPtr p_neighbour_cell_i = rCellPopulation.GetCellUsingLocationIndex(*neighbouring_node_indices.begin());
        unsigned second_neighbour_indice = GetNthNeighbourIndice(neighbouring_node_indices, 1);
        CellPtr p_neighbour_cell_k = rCellPopulation.GetCellUsingLocationIndex(second_neighbour_indice);
        c_vector<double,DIM> xi = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_i);
        c_vector<double,DIM> xk = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_k);
    
        u = rCellPopulation.rGetMesh().GetVectorFromAtoB(xi, xj);
        v = rCellPopulation.rGetMesh().GetVectorFromAtoB(xk, xj);

        // alphangularmin needs to be in between 0 and pi 
        double alphangularmin = GetAngleFromVectors(u, v);
        alphangularmin = fmod(alphangularmin, 2*M_PI);
        if(alphangularmin < 0){
            alphangularmin += 2*M_PI;
        }
        if(alphangularmin > M_PI){
            alphangularmin = 2*M_PI - alphangularmin;
        }

        alphangular = alphangularmin;

        if(neighbouring_node_indices.size() > 2){
            double alpha;
            // if there are more than two neighbours, then we need to sort out which one are making the smallest angle 
            c_vector<double, DIM> ximin = xi;
            c_vector<double, DIM> xkmin = xk;
            
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
                        
                        // alpha needs to be between 0 and pi
                        alpha = GetAngleFromVectors(u, v);
                        alpha = fmod(alpha, 2*M_PI);
                        if(alpha < 0){
                            alpha += 2*M_PI;
                        }
                        if(alpha > M_PI){
                            alpha = 2*M_PI - alpha;
                        }

                        if(alpha < alphangularmin){
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
        }
    }

    return std::make_tuple(alphangular, u, v);
}

// overrides AddForceContribution()
template<unsigned DIM>
void AngularForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{    
    NodeBasedCellPopulation<DIM>* p_node_population = dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);
    if (!p_node_population)
    {
        EXCEPTION("NodeBasedCellPopulation only valid for node.");
    }

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin(); cell_iter != rCellPopulation.End(); ++cell_iter)
    {
        if (cell_iter->GetMutationState()->template IsType<VesselSegmentMutationState>())
        {
            c_vector<double, DIM> angularforce = zero_vector<double>(DIM);

            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            std::set<unsigned> neighbouring_node_indices = GetSameBranchNeighbours(rCellPopulation,p_node_population, *cell_iter);

            // we calculate the two neigbours and the angle they make with the node considered
            double alphangular;
            c_vector<double, DIM> u;
            c_vector<double, DIM> v;
            //std::tie(alphangular, u, v) = ClosestAngleVesselSegment(rCellPopulation, (*cell_iter), neighbouring_node_indices);
            std::tie(alphangular, u, v) = OptimalAngleVesselElement(rCellPopulation, (*cell_iter), neighbouring_node_indices);

            c_vector<double, DIM> r_angularneighbours = u + v;
            double magnitude_of_angularneighbours = norm_2(r_angularneighbours);

            if(magnitude_of_angularneighbours != 0){
                angularforce = mOmegaa*(M_PI-alphangular)/magnitude_of_angularneighbours*r_angularneighbours;
            }

            rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(angularforce);
            
            // test time force
            // for(int i = 0; i < 100; ++i){rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(angularforce);}
        }
    }
}

template<unsigned DIM>
void AngularForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<Angular Persistence Spring Constant>" << mOmegaa << "</Omegaa>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class AngularForce<1>;
template class AngularForce<2>;
template class AngularForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AngularForce)