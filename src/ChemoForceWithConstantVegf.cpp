#include "ChemoForceWithConstantVegf.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "VegfEquationPde.hpp"
#include "LinearBasisFunction.hpp"
#include "ReplicatableVector.hpp"

#include "Debug.hpp"

template<unsigned DIM>
ChemoForceWithConstantVegf<DIM>::ChemoForceWithConstantVegf(double chiConstantVegf, double cx)
    : ChemoForce<DIM>(chiConstantVegf, 0.0, 0.0, 0.0), mChiConstantVegf(chiConstantVegf), mCX(cx)
{
    assert(chiConstantVegf>0);
}

template<unsigned DIM>
ChemoForceWithConstantVegf<DIM>::~ChemoForceWithConstantVegf()
{
}

template<unsigned DIM>
c_vector<double, DIM>& ChemoForceWithConstantVegf<DIM>::GetGradient(unsigned node_index)
{
    return mGradientsVegfAnalyticalApproxPde[node_index];
}

template<unsigned DIM>
double ChemoForceWithConstantVegf<DIM>::GetMagnitudeGradient(unsigned node_index)
{
    return norm_2(GetGradient(node_index))/mChiAnalyticalApproxPde;
}

template<unsigned DIM>
void ChemoForceWithConstantVegf<DIM>::CalculateVegfGradient(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Initialise gradients size
    unsigned num_nodes = rCellPopulation.GetNumNodes();
    mGradientsVegfAnalyticalApproxPde.resize(num_nodes, zero_vector<double>(DIM));

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin(); cell_iter != rCellPopulation.End(); ++cell_iter)
    {
        // we collect the cell data necessary (node index and cell pointer)
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        CellPtr pCell = rCellPopulation.GetCellUsingLocationIndex(node_index); 

        c_vector<double, DIM> x_parent = rCellPopulation.GetLocationOfCellCentre(pCell);

        c_vector<double, DIM> r_gradient_cell = zero_vector<double>(DIM);

        if (pCell->GetMutationState()->IsType<TipCellMutationState>())
        {
            if(DIM == 3){
                r_gradient_cell(0) = -mCX; 
                r_gradient_cell(1) = 0.0; 
                r_gradient_cell(2) = 0.0; 
            } else if (DIM == 2){
                r_gradient_cell(0) = -mCX; 
                r_gradient_cell(1) = 0.0; 
            } else {
                r_gradient_cell(0) = -mCX; 
            }

            mGradientsVegfAnalyticalApproxPde[node_index] += r_gradient_cell;
        }
    }
}

// Explicit instantiation
template class ChemoForceWithConstantVegf<1>;
template class ChemoForceWithConstantVegf<2>;
template class ChemoForceWithConstantVegf<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemoForceWithConstantVegf)