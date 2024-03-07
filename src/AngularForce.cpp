#include "AngularForce.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "NodeBasedCellPopulation.hpp"

#include <tuple>
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
std::tuple<double, c_vector<double,DIM>, c_vector<double,DIM>> AngularForce<DIM>::GetAngleVesselElement(CellPtr cell_ptr, std::set<unsigned> neighbouring_node_indices, AbstractCellPopulation<DIM>& rCellPopulation)
{
    c_vector<double,DIM> xj = rCellPopulation.GetLocationOfCellCentre(cell_ptr);

    double alphangular = M_PI; // we initialise such as the force is equal to 0, geometrically, it means that the elements are aligned
    c_vector<double, DIM> u;
    c_vector<double, DIM> v;

    CellPtr p_neighbour_cell_i = rCellPopulation.GetCellUsingLocationIndex(*neighbouring_node_indices.begin());
    CellPtr p_neighbour_cell_k = rCellPopulation.GetCellUsingLocationIndex(*neighbouring_node_indices.begin()+1);
    c_vector<double,DIM> xi = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_i);
    c_vector<double,DIM> xk = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_k);
    
    u = rCellPopulation.rGetMesh().GetVectorFromAtoB(xi, xj);
    v = rCellPopulation.rGetMesh().GetVectorFromAtoB(xk, xj);

    double scalar_product_uv;
    if(DIM == 2){
        scalar_product_uv = u[0]*v[0] + u[1]*v[1];
    } else if(DIM == 3){
        scalar_product_uv = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
    } else{
        scalar_product_uv = u[0]*v[0];
    }

    double alphangularmin = 0;
    if(norm_2(u)!= 0 && norm_2(v) != 0){
        alphangularmin = std::acos(scalar_product_uv/(norm_2(u)*norm_2(v)));
    } 

    if(neighbouring_node_indices.size() > 2){
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
                    double scalar_product_xijk;
                    if(DIM == 2){
                        scalar_product_xijk = u[0]*v[0] + u[1]*v[1];
                    } else if(DIM == 3){
                        scalar_product_xijk = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
                    } else{
                        scalar_product_xijk = u[0]*v[0];
                    }                        
                    double alpha = 0;
                    if(norm_2(xij) != 0 && norm_2(xkj) != 0){
                        alpha = std::acos(scalar_product_xijk/(norm_2(xij)*norm_2(xkj)));
                    }

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
        double scalar_product_uv;
        if(DIM == 2){
            scalar_product_uv = u[0]*v[0] + u[1]*v[1];
        } else if(DIM == 3){
            scalar_product_uv = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
        } else{
            scalar_product_uv = u[0]*v[0];
        }
        if(norm_2(u)!=0 && norm_2(v)!=0){
            alphangular = std::acos(scalar_product_uv/(norm_2(u)*norm_2(v)));
        }
    } else if(neighbouring_node_indices.size() == 2){
        alphangular = alphangularmin;
    } 

    return std::make_tuple(alphangular, u, v);
}

// overrides AddForceContribution()
template<unsigned DIM>
void AngularForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    //TRACE("Begin Angular Force");
    // we create the vector force 
    c_vector<double, DIM> angularforce = zero_vector<double>(DIM); 
    
    NodeBasedCellPopulation<DIM>* p_node_population = dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);
    if (!p_node_population)
    {
        EXCEPTION("NodeBasedCellPopulation only valid for node.");
    }

    for (unsigned node_index=0; node_index<rCellPopulation.GetNumNodes(); node_index++)
    {
        CellPtr cell_ptr = rCellPopulation.GetCellUsingLocationIndex(node_index);
        if (cell_ptr->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
        {
            std::set<unsigned> neighbouring_node_indices = p_node_population->GetNeighbouringNodeIndices(node_index);

            //PRINT_VARIABLE(neighbouring_node_indices.size());

            if(neighbouring_node_indices.size() >= 2){
                double alphangular;
                c_vector<double, DIM> u;
                c_vector<double, DIM> v;
                std::tie(alphangular, u, v) = GetAngleVesselElement(cell_ptr, neighbouring_node_indices, rCellPopulation);
                //PRINT_VARIABLE(alphangular);

                c_vector<double, DIM> r_angularneighbours = u + v;
                double magnitude_of_angularneighbours = norm_2(r_angularneighbours);

                c_vector<double,DIM> angularforce;
                if(magnitude_of_angularneighbours != 0){
                    angularforce = -mOmegaa*(M_PI-alphangular)/magnitude_of_angularneighbours*r_angularneighbours;
                } else{
                    angularforce = zero_vector<double>(DIM);
                }
                rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(angularforce);
            } else {
                c_vector<double,DIM> angularforce = zero_vector<double>(DIM);
                rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(angularforce);
            }
        }
    }
    //TRACE("End");
}

template<unsigned DIM>
void AngularForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    double omegaa = GetAngularPersistence();
    *rParamsFile << "\t\t\t<Angular Persistence Spring Constant>" << omegaa << "</Omegaa>\n";

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
//CHASTE_CLASS_EXPORT(AngularForce)