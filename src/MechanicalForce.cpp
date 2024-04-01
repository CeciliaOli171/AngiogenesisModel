#include "MechanicalForce.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "NodeBasedCellPopulation.hpp"

template<unsigned DIM>
MechanicalForce<DIM>::MechanicalForce(double Sc)
    : AbstractForce<DIM>()
{
    assert(Sc > 0);
    mSc = Sc;
}

template<unsigned DIM>
MechanicalForce<DIM>::~MechanicalForce()
{
}

template<unsigned DIM>
double MechanicalForce<DIM>::GetSpringConstant()
{ 
    return mSc;
}

// function that returns the approximation of the length of the vessel element k
template<unsigned DIM>
double MechanicalForce<DIM>::LengthVesselSegment(AbstractCellPopulation<DIM, DIM>& rCellPopulation, CellPtr pCell, std::set<unsigned> neighbouring_node_indices){
    // initialisation
    double length_vessel_segment = 0; 

    // we select the coordinates of the element 
    c_vector<double,DIM> xj = rCellPopulation.GetLocationOfCellCentre(pCell);

    // we consider the neighbours of the element 
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

// Spring function depending on the distance between the cells 
template<unsigned DIM>
double MechanicalForce<DIM>::SpringFunction(double x, double lc, double Rc){
    double S;

    if(x > 0){
        S = mSc*x; // compressed 
    } else if (x <= 0 && x > 2*Rc-lc){ 
        S = 100*mSc*x; // elongated 
    } else {
        S = 0;
    }

    return S;
}

// overrides AddForceContribution()
template<unsigned DIM>
void MechanicalForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // initialisation 
    c_vector<double, DIM> mechanicalforce;

    NodeBasedCellPopulation<DIM>* p_node_population = dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);
    if (!p_node_population)
    {
        EXCEPTION("NodeBasedCellPopulation only valid for node");
    }
 
    // we applied the force to the cell population corresponding 
    // here, the random force is applied to every cell 
    unsigned node_index = 0;
    for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
         node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        // we collect cell data (pointer, indice, coordinates, neighbours, radius, length)
        CellPtr pCell = rCellPopulation.GetCellUsingLocationIndex(node_index);
        Node<DIM>* pNodeCell = rCellPopulation.GetNode(node_index);
        c_vector<double,DIM> xk = rCellPopulation.GetLocationOfCellCentre(pCell);
        std::set<unsigned> neighbouring_node_indices = p_node_population->GetNeighbouringNodeIndices(node_index);
        double Rc = pNodeCell->GetRadius();
        double lc = LengthVesselSegment(rCellPopulation, pCell, neighbouring_node_indices);

        mechanicalforce = zero_vector<double>(DIM);

        for(std::set<unsigned>::iterator j = neighbouring_node_indices.begin();
        j != neighbouring_node_indices.end();
        ++j){
            CellPtr p_neighbour_cell_j = rCellPopulation.GetCellUsingLocationIndex(*j);
            c_vector<double, DIM> xj = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell_j);
            c_vector<double, DIM> u = xk - xj;
            double magnitude_u = norm_2(u);
            double S = SpringFunction(2*Rc-magnitude_u, lc, Rc);
            mechanicalforce = mechanicalforce + S*u/magnitude_u;
        }

        node_iter->AddAppliedForceContribution(mechanicalforce);
        ++node_index;
    }
}

template<unsigned DIM>
void MechanicalForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<Spring Constant>" << mSc << "</Sc>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class MechanicalForce<1>;
template class MechanicalForce<2>;
template class MechanicalForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MechanicalForce)