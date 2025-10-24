#include "ChemoForceWithPdes.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "VegfEquationPde.hpp"
#include "LinearBasisFunction.hpp"
#include "ReplicatableVector.hpp"

#include "Debug.hpp"

template<unsigned DIM>
ChemoForceWithPdes<DIM>::ChemoForceWithPdes(double chiPdes, double hx, double hy, double hz, boost::shared_ptr<AbstractPdeModifier<DIM> > pPdeModifier)
    : ChemoForce<DIM>(chiPdes, hx, hy, hz, 0.0, 0.0, 0.0), mChiPdes(chiPdes), mHX(hx), mHY(hy), mHZ(hz),  mpPdeModifier(pPdeModifier)
{
    assert(chiPdes>0);
}

template<unsigned DIM>
ChemoForceWithPdes<DIM>::~ChemoForceWithPdes()
{
}

template<unsigned DIM>
c_vector<double, DIM>& ChemoForceWithPdes<DIM>::GetGradient(AbstractCellPopulation<DIM>& rCellPopulation, unsigned node_index)
{
    TetrahedralMesh<DIM,DIM>* r_mesh = mpPdeModifier->GetFeMesh();
    ReplicatableVector previous_solution(mpPdeModifier->GetSolution());

    CellPtr pCell = rCellPopulation.GetCellUsingLocationIndex(node_index);
    unsigned elem_index = mpPdeModifier->GetFeMesh()->GetContainingElementIndex(rCellPopulation.GetLocationOfCellCentre(pCell));
    Element<DIM,DIM>& r_elem = *(r_mesh->GetElement(elem_index)); // element where the vessel tip is located 

    // Calculate the basis functions at any point (eg zero) in the element
    c_matrix<double, DIM, DIM> jacobian, inverse_jacobian;
    double jacobian_det;
    const ChastePoint<DIM> zero_point;
    c_matrix<double, DIM, DIM+1> grad_phi;
    c_vector<double, DIM> gradient_on_element;

    r_mesh->GetInverseJacobianForElement(elem_index, jacobian, jacobian_det, inverse_jacobian);
    LinearBasisFunction<DIM>::ComputeTransformedBasisFunctionDerivatives(zero_point, inverse_jacobian, grad_phi);

    for (unsigned node_index=0; node_index<DIM+1; node_index++)
    {
        unsigned node_global_index = r_elem.GetNodeGlobalIndex(node_index);
        double pde_solution = previous_solution[node_global_index]; // Get PDE solution

        // Interpolate gradient
        for (unsigned i=0; i<DIM; i++)
        {
            gradient_on_element(i) += pde_solution* grad_phi(i, node_index);
        }
    }

    // add chemoattractant 
    c_vector<double, DIM> r_gradient_cell = gradient_on_element;
    if(DIM == 3){
        r_gradient_cell(0) = mChiPdes*(r_gradient_cell(0)-mHX); 
        r_gradient_cell(1) = mChiPdes*(r_gradient_cell(1)-mHY); 
        r_gradient_cell(2) = mChiPdes*(r_gradient_cell(2)-mHZ); 
    } else if (DIM == 2){
        r_gradient_cell(0) = mChiPdes*(r_gradient_cell(0)-mHX); 
        r_gradient_cell(1) = mChiPdes*(r_gradient_cell(1)-mHY); 
    } else {
        r_gradient_cell(0) = mChiPdes*(r_gradient_cell(0)-mHX); 
    }

    return r_gradient_cell;
}

// Explicit instantiation
template class ChemoForceWithPdes<1>;
template class ChemoForceWithPdes<2>;
template class ChemoForceWithPdes<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemoForceWithPdes)