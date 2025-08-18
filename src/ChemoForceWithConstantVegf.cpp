#include "ChemoForceWithConstantVegf.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "VegfEquationPde.hpp"
#include "LinearBasisFunction.hpp"
#include "ReplicatableVector.hpp"

#include "Debug.hpp"

template<unsigned DIM>
ChemoForceWithConstantVegf<DIM>::ChemoForceWithConstantVegf(double chiConstantVegf, double hx, double constantBackground)
    : ChemoForce<DIM>(chiConstantVegf, hx, 0.0, 0.0, 0.0, 0.0, 0.0), mChiConstantVegf(chiConstantVegf), mConstantBackground(constantBackground)
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
    return norm_2(GetGradient(node_index));
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
        if (cell_iter->GetMutationState()->template IsType<TipCellMutationState>())
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

            c_vector<double, DIM> r_gradient_cell = zero_vector<double>(DIM);
            
            r_gradient_cell(0) = 0.0; 

            mGradientsVegfAnalyticalApproxPde[node_index] = r_gradient_cell;
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