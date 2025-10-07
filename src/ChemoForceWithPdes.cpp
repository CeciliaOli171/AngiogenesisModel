#include "ChemoForceWithPdes.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "VegfEquationPde.hpp"
#include "LinearBasisFunction.hpp"
#include "ReplicatableVector.hpp"

#include "Debug.hpp"

template<unsigned DIM>
ChemoForceWithPdes<DIM>::ChemoForceWithPdes(double chiPdes, double hx, double hy, double hz, boost::shared_ptr<AbstractBoxDomainPdeModifier<DIM> > pPdeModifier)
    : ChemoForce<DIM>(chiPdes, hx, hy, hz, 0.0, 0.0, 0.0), mChiPdes(chiPdes),  mpPdeModifier(pPdeModifier)
{
    assert(chiPdes>0);
}

template<unsigned DIM>
ChemoForceWithPdes<DIM>::~ChemoForceWithPdes()
{
}

template<unsigned DIM>
c_vector<double, DIM>& ChemoForceWithPdes<DIM>::GetGradient(unsigned node_index)
{
    return mGradientsVegfPdes[node_index];
}

template<unsigned DIM>
double ChemoForceWithPdes<DIM>::GetMagnitudeGradient(unsigned node_index)
{
    return norm_2(GetGradient(node_index));
}

template<unsigned DIM>
void ChemoForceWithPdes<DIM>::CalculateVegfGradient(AbstractCellPopulation<DIM>& rCellPopulation)
{
    TetrahedralMesh<DIM,DIM>* r_mesh = mpPdeModifier->GetFeMesh();
    ReplicatableVector previous_solution(mpPdeModifier->GetSolution());
 
    // Initialise gradients size
    unsigned num_nodes = rCellPopulation.GetNumNodes();
    mGradientsVegfPdes.resize(num_nodes, zero_vector<double>(DIM));

    // The constant gradients at each element
    std::vector<c_vector<double, DIM> > gradients_on_elements;
    unsigned num_elements = r_mesh->GetNumElements();
    gradients_on_elements.resize(num_elements, zero_vector<double>(DIM));

    for (unsigned elem_index=0; elem_index<num_elements; elem_index++)
    {
        Element<DIM,DIM>& r_elem = *(r_mesh->GetElement(elem_index));

        // Calculate the basis functions at any point (eg zero) in the element
        c_matrix<double, DIM, DIM> jacobian, inverse_jacobian;
        double jacobian_det;
        r_mesh->GetInverseJacobianForElement(elem_index, jacobian, jacobian_det, inverse_jacobian);
        const ChastePoint<DIM> zero_point;
        c_matrix<double, DIM, DIM+1> grad_phi;
        LinearBasisFunction<DIM>::ComputeTransformedBasisFunctionDerivatives(zero_point, inverse_jacobian, grad_phi);

        for (unsigned node_index=0; node_index<DIM+1; node_index++)
        {
            unsigned node_global_index = r_elem.GetNodeGlobalIndex(node_index);

            // Get PDE solution
            double pde_solution = previous_solution[node_global_index];

            // Interpolate gradient
            for (unsigned i=0; i<DIM; i++)
            {
                gradients_on_elements[elem_index](i) += pde_solution* grad_phi(i, node_index);
            }
        }
    }

    // we need to go through all the cells and collect the value of the gradients at those location
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        unsigned node_indice = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

        unsigned elem_index = mpPdeModifier->GetFeMesh()->GetContainingElementIndex(rCellPopulation.GetLocationOfCellCentre((*cell_iter)));

        mGradientsVegfPdes[node_indice] += mChiPdes*gradients_on_elements[elem_index];
    }
}

// Explicit instantiation
template class ChemoForceWithPdes<1>;
template class ChemoForceWithPdes<2>;
template class ChemoForceWithPdes<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemoForceWithPdes)