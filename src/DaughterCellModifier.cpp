#include "DaughterCellModifier.hpp"
#include "MeshBasedCellPopulation.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"

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
double DaughterCellModifier<DIM>::OptimalAngleVesselSegment(std::set<unsigned> neighbouring_node_indices, CellPtr cell_ptr, AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    c_vector<double, DIM> xj = rCellPopulation.GetLocationOfCellCentre(cell_ptr);

    double alphangular = M_PI; // we initialise such as the force is equal to 0 
    c_vector<double, DIM> u;
    c_vector<double, DIM> v;

    CellPtr p_neighbour_cell_i = rCellPopulation.GetCellUsingLocationIndex(*neighbouring_node_indices.begin());
    CellPtr p_neighbour_cell_k = rCellPopulation.GetCellUsingLocationIndex(*neighbouring_node_indices.end());
    c_vector<double,DIM> xi = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_i);
    c_vector<double,DIM> xk = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_k);
    
    u = rCellPopulation.rGetMesh().GetVectorFromAtoB(xi, xj);
    v = rCellPopulation.rGetMesh().GetVectorFromAtoB(xk, xj);
    double scalar_product_uv;
    double alphangularmin = 0;
    if(DIM == 2){
        scalar_product_uv = u[0]*v[0] + u[1]*v[1];
    } else if(DIM == 3){
        scalar_product_uv = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
    } else{
        scalar_product_uv = u[0]*v[0];
    }
    if(norm_2(u) !=0 && norm_2(v) != 0){
        alphangularmin = std::acos(scalar_product_uv/(norm_2(u)*norm_2(v)));
    } else {
        TRACE("0 norm found");
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
        double scalar_product_uv;
        if(DIM == 2){
            scalar_product_uv = u[0]*v[0] + u[1]*v[1];
        } else if(DIM == 3){
            scalar_product_uv = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
        } else{
            scalar_product_uv = u[0]*v[0];
        }
        if(norm_2(u) != 0 && norm_2(v) != 0){
            alphangular = std::acos(scalar_product_uv/(norm_2(u)*norm_2(v)));
        }
        cout << alphangular << endl; 
    } else if(neighbouring_node_indices.size() == 2){
        alphangular = alphangularmin;
    } 
    PRINT_VARIABLE(alphangular);
    return alphangular;
}

// function that returns the approximation of the length of the vessel element k
template<unsigned DIM>
double DaughterCellModifier<DIM>::LengthVesselSegment(std::set<unsigned> neighbouring_node_indices, CellPtr cell_ptr, AbstractCellPopulation<DIM, DIM>& rCellPopulation){
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

    // cout << "Size of new cell population = " << rCellPopulation.GetNumNodes() << endl;
    NodeBasedCellPopulation<DIM>* p_node_population = dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);
    if (!p_node_population)
    {
        EXCEPTION("NodeBasedCellPopulation only valid for node");
    }

    // we go through all the cells from the new cell population
    for(unsigned node_index = mOldNumNodes; node_index<rCellPopulation.GetNumNodes(); node_index++){
        // if the old cell population does not have the index then we have the daughter cell's position 
        CellPtr pDaughterCell = rCellPopulation.GetCellUsingLocationIndex(node_index);

        // we check the angle made by the daughter cell with its neighbours
        unsigned daughter_node_index = node_index; // node index
        Node<DIM>* p_daughter_cell = rCellPopulation.GetNode(daughter_node_index);
        std::set<unsigned> neighbouring_node_indices = p_node_population->GetNeighbouringNodeIndices(daughter_node_index);

        double angle_vessel_segment = OptimalAngleVesselSegment(neighbouring_node_indices, pDaughterCell, rCellPopulation);
        double random_angle = 2.0*M_PI*RandomNumberGenerator::Instance()->ranf(); 
        cout << "random angle = " << random_angle << endl;

        // depending on the angle made by the random vector and the parent : 
        // the new cell will start a new sprout and differentiate into a tip cell 
        // or it will be incorporated into the parent vessel and differentiate into a vessel segment 

        // sprouting probability : to change because all the new cells are currently vessel segment 
        // double length_vessel_segment = LengthVesselSegment(neighbouring_node_indices, pDaughterCell, rCellPopulation);
        // double Rc = p_daughter_cell->GetRadius();
        // double kspr = 2e-4;
        // double Kspr = 0.01;
        // double Psprout = kspr/(Kspr + std::max(2*Rc-length_vessel_segment, 0.0));
        
        double Psprout = 0.5;
        double rand_number = RandomNumberGenerator::Instance()->ranf();
        // cout << "random number =" << rand_number << endl;

        if(rand_number < Psprout){ //random_angle >= (1-Psprout)*angle_vessel_segment/2 && random_angle <= (1+Psprout)*angle_vessel_segment/2){
            // the new cell is a tip cell 
            boost::shared_ptr<AbstractCellProperty> p_tipcell_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
            pDaughterCell->SetCellProliferativeType(p_tipcell_type);
            cout << "test : new cell is a tip cell" << endl;
        } else {
            // the new cell is a vessel segment 
            boost::shared_ptr<AbstractCellProperty> p_vesselsegment_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
            pDaughterCell->SetCellProliferativeType(p_vesselsegment_type);
            cout << "test : new cell is a vessel segment" << endl;
        }
    }

    mOldNumNodes = rCellPopulation.GetNumNodes();

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

