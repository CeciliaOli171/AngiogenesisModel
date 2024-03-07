#include "Sprouting.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"

#include "RandomNumberGenerator.hpp"
#include <algorithm>
#include "Debug.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SproutingRule<ELEMENT_DIM, SPACE_DIM>::SproutingRule()
    : AbstractCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SproutingRule<ELEMENT_DIM, SPACE_DIM>::~SproutingRule()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double SproutingRule<ELEMENT_DIM, SPACE_DIM>::OptimalAngleVesselSegment(std::set<unsigned> neighbouring_node_indices, CellPtr cell_ptr, AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    c_vector<double,SPACE_DIM> xj = rCellPopulation.GetLocationOfCellCentre(cell_ptr);

    double alphangular = M_PI; // we initialise such as the force is equal to 0 
    c_vector<double, SPACE_DIM> u;
    c_vector<double, SPACE_DIM> v;

    CellPtr p_neighbour_cell_i = rCellPopulation.GetCellUsingLocationIndex(*neighbouring_node_indices.begin());
    CellPtr p_neighbour_cell_k = rCellPopulation.GetCellUsingLocationIndex(*neighbouring_node_indices.end());
    c_vector<double,SPACE_DIM> xi = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_i);
    c_vector<double,SPACE_DIM> xk = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_k);
    
    u = rCellPopulation.rGetMesh().GetVectorFromAtoB(xi, xj);
    v = rCellPopulation.rGetMesh().GetVectorFromAtoB(xk, xj);
    double scalar_product_uv;
    double alphangularmin = 0;
    if(SPACE_DIM == 2){
        scalar_product_uv = u[0]*v[0] + u[1]*v[1];
    } else if(SPACE_DIM == 3){
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
        c_vector<double, SPACE_DIM> ximin = xi;
        c_vector<double, SPACE_DIM> xkmin = xk;
        
        for(std::set<unsigned>::iterator i = neighbouring_node_indices.begin();
        i != neighbouring_node_indices.end();
        ++i){
            for(std::set<unsigned>::iterator k = neighbouring_node_indices.begin();
            k != neighbouring_node_indices.end();
            ++k){
                if(i != k){
                    CellPtr p_neighbour_cell_i = rCellPopulation.GetCellUsingLocationIndex(*i);
                    CellPtr p_neighbour_cell_k = rCellPopulation.GetCellUsingLocationIndex(*k);
                    c_vector<double, SPACE_DIM> xi = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_i);
                    c_vector<double, SPACE_DIM> xk = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_k);
                    c_vector<double, SPACE_DIM> xij = rCellPopulation.rGetMesh().GetVectorFromAtoB(xi, xj);
                    c_vector<double, SPACE_DIM> xkj = rCellPopulation.rGetMesh().GetVectorFromAtoB(xk, xj);
                    double scalar_product_xijk;
                    if(SPACE_DIM == 2){
                        scalar_product_xijk = u[0]*v[0] + u[1]*v[1];
                    } else if(SPACE_DIM == 3){
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
        if(SPACE_DIM == 2){
            scalar_product_uv = u[0]*v[0] + u[1]*v[1];
        } else if(SPACE_DIM == 3){
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
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double SproutingRule<ELEMENT_DIM, SPACE_DIM>::LengthVesselSegment(std::set<unsigned> neighbouring_node_indices, CellPtr cell_ptr, AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation){
    double length_vessel_segment = 0; 

    // select the element 
    c_vector<double,SPACE_DIM> xj = rCellPopulation.GetLocationOfCellCentre(cell_ptr);

    // consider the neighbours of the element 
    if(neighbouring_node_indices.size() > 2){ 
        for(std::set<unsigned>::iterator k = neighbouring_node_indices.begin();
        k != neighbouring_node_indices.end();
        ++k){
            CellPtr p_neighbour_cell_k = rCellPopulation.GetCellUsingLocationIndex(*k);
            c_vector<double, SPACE_DIM> xk = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_k);

            c_vector<double, SPACE_DIM> u = rCellPopulation.rGetMesh().GetVectorFromAtoB(xk, xj);
            length_vessel_segment = length_vessel_segment + norm_2(u);
        }
        length_vessel_segment = length_vessel_segment/neighbouring_node_indices.size();
    } else if(neighbouring_node_indices.size() == 2){
        CellPtr p_neighbour_cell_i = rCellPopulation.GetCellUsingLocationIndex(*neighbouring_node_indices.begin());
        CellPtr p_neighbour_cell_k = rCellPopulation.GetCellUsingLocationIndex(*neighbouring_node_indices.end());
        c_vector<double,SPACE_DIM> xi = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_i);
        c_vector<double,SPACE_DIM> xk = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_k);

        c_vector<double, SPACE_DIM> u = rCellPopulation.rGetMesh().GetVectorFromAtoB(xi, xj);
        c_vector<double, SPACE_DIM> v = rCellPopulation.rGetMesh().GetVectorFromAtoB(xk, xj);
        length_vessel_segment = (norm_2(u)+norm_2(v))/2; 
    }

    // depending on the number of neighbours, we calculate the length of the vessel element
    return length_vessel_segment;
}

// function that set the new daughter cell property : tip cell or vessel segment
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SproutingRule<ELEMENT_DIM, SPACE_DIM>::DaughterTypeOfCell(AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, unsigned OldNumNodes){

    cout << "Size of new cell population = " << rCellPopulation.GetNumNodes() << endl;
    NodeBasedCellPopulation<SPACE_DIM>* p_node_population = dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation);
    if (!p_node_population)
    {
        EXCEPTION("NodeBasedCellPopulation only valid for node");
    }

    // we go through all the cells from the new cell population
    for(unsigned node_index = OldNumNodes; node_index<rCellPopulation.GetNumNodes(); node_index++){
        // if the old cell population does not have the index then we have the daughter cell's position 
        CellPtr pDaughterCell = rCellPopulation.GetCellUsingLocationIndex(node_index);

        // we check the angle made by the daughter cell with its neighbours
        unsigned daughter_node_index = node_index; // node index
        Node<SPACE_DIM>* p_daughter_cell = rCellPopulation.GetNode(daughter_node_index);
        std::set<unsigned> neighbouring_node_indices = p_node_population->GetNeighbouringNodeIndices(daughter_node_index);

        double angle_vessel_segment = OptimalAngleVesselSegment(neighbouring_node_indices, pDaughterCell, rCellPopulation);
        double random_angle = 2.0*M_PI*RandomNumberGenerator::Instance()->ranf(); 
        cout << "random angle = " << random_angle << endl;

        // depending on the angle made by the random vector and the parent : 
        // the new cell will start a new sprout and differentiate into a tip cell 
        // or it will be incorporated into the parent vessel and differentiate into a vessel segment 

        // sprouting probability : to change because all the new cells are currently vessel segment 
        double length_vessel_segment = LengthVesselSegment(neighbouring_node_indices, pDaughterCell, rCellPopulation);
        double Rc = p_daughter_cell->GetRadius();
        double kspr = 2e-4;
        double Kspr = 0.01;
        double Psprout = kspr/(Kspr + std::max(2*Rc-length_vessel_segment, 0.0));
        
        Psprout = 0.5;
        double rand_number = RandomNumberGenerator::Instance()->ranf();
        cout << "random number =" << rand_number << endl;

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

}

// Overridden CalculateCellDivisionVector() method, return the daughter cell positions
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > SproutingRule<ELEMENT_DIM, SPACE_DIM>::CalculateCellDivisionVector(CellPtr pParentCell, AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    TRACE("Begin");
    // consider a parent cell 
    c_vector<double, SPACE_DIM> parent_position = rCellPopulation.GetLocationOfCellCentre(pParentCell);

    // make a random direction vector of the required length for the daughter cell 
    c_vector<double, SPACE_DIM> random_vector; 

    double random_angle = 2.0*M_PI*RandomNumberGenerator::Instance()->ranf();

    if(SPACE_DIM == 2){
        random_vector(0) = 0.5*cos(random_angle);
        random_vector(1) = 0.5*sin(random_angle);
    } else if(SPACE_DIM == 3){
        double random_angle_2 = 2.0*M_PI*RandomNumberGenerator::Instance()->ranf();
        random_vector(0) = 0.5*sin(random_angle_2)*cos(random_angle);
        random_vector(1) = 0.5*sin(random_angle_2)*sin(random_angle);
        random_vector(2) = 0.5*cos(random_angle);
    } else{
        random_vector(0) = 0.5*random_angle;
    }

    NodeBasedCellPopulation<SPACE_DIM>* p_node_population = dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation);
    if (!p_node_population)
    {
        EXCEPTION("NodeBasedCellPopulation only valid for node");
    }

    unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(pParentCell); // node indice 
    std::set<unsigned> neighbouring_node_indices = p_node_population->GetNeighbouringNodeIndices(node_index);

    // we consider the length of the vessel segment and the angle between the vessel segment and its two closest neighbours 
    //double length_vessel_segment = LengthVesselSegment(neighbouring_node_indices, pParentCell, rCellPopulation);
    //double angle_vessel_segment = OptimalAngleVesselSegment(neighbouring_node_indices, pParentCell, rCellPopulation); 

    // creation of the daughter cell to the map
    c_vector<double, SPACE_DIM> daughter_position = parent_position + 3.0*random_vector; 

    if(SPACE_DIM == 2){
        cout << "parent position = " << parent_position(0) << ';' << parent_position(1) << endl;
        cout << "daughter position = "<< daughter_position(0) << ';' << daughter_position(1) << endl;
        cout << " " << endl;
    } else if(SPACE_DIM == 3){
        cout << "parent position = " << parent_position(0) << ';' << parent_position(1) << ';' << parent_position(2) << endl;
        cout << "daughter position = " << daughter_position(0) << ';' << daughter_position(1) << ';' << daughter_position(2) << endl;
        cout << " " << endl;
    } else{
        cout << "parent position = " << parent_position(0) << endl;
        cout << "daughter position = " << daughter_position(0) << endl;
        cout << " " << endl;
    }

    std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > positions(parent_position, daughter_position);

    TRACE("End");

    return positions;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SproutingRule<ELEMENT_DIM, SPACE_DIM>::OutputCellCentreBasedDivisionRuleInfo(out_stream& rParamsFile){
    std::string cell_division_rule_type = "Sprouting Division Rule";

    *rParamsFile << "\t\t\t<" << cell_division_rule_type << ">\n";

    // Do we need it ?
    //OutputCellCentreBasedDivisionRuleParameters(rParamsFile);
    //*rParamsFile << "\t\t\t</" << cell_division_rule_type << ">\n";
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