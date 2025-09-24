#include "ChemoForceWithAnalyticalApproximationPde.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "VegfEquationPde.hpp"
#include "LinearBasisFunction.hpp"
#include "ReplicatableVector.hpp"

#include "Debug.hpp"

template<unsigned DIM>
ChemoForceWithAnalyticalApproximationPde<DIM>::ChemoForceWithAnalyticalApproximationPde(double chiAnalyticalApproxPde, double hx, double diffusionCoefficient, double decayCoefficient, double creationCoefficient, double consumptionCoefficient, double sourceValue, double constantBackground)
    : ChemoForce<DIM>(chiAnalyticalApproxPde, hx, 0.0, 0.0, 0.0, 0.0, 0.0), mChiAnalyticalApproxPde(chiAnalyticalApproxPde), mDiffusionCoefficient(diffusionCoefficient), mDecayCoefficient(decayCoefficient), mCreationCoefficient(creationCoefficient), mConsumptionCoefficient(consumptionCoefficient), mSourceValue(sourceValue), mConstantBackground(constantBackground)
{
    assert(chiAnalyticalApproxPde>0);
}

template<unsigned DIM>
ChemoForceWithAnalyticalApproximationPde<DIM>::~ChemoForceWithAnalyticalApproximationPde()
{
}

template<unsigned DIM>
c_vector<double, DIM>& ChemoForceWithAnalyticalApproximationPde<DIM>::GetGradient(unsigned node_index)
{
    return mGradientsVegfAnalyticalApproxPde[node_index];
}

template<unsigned DIM>
double ChemoForceWithAnalyticalApproximationPde<DIM>::GetMagnitudeGradient(unsigned node_index)
{
    return norm_2(GetGradient(node_index));
}

template<unsigned DIM>
void ChemoForceWithAnalyticalApproximationPde<DIM>::CalculateVegfGradient(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Initialise gradients size
    unsigned num_nodes = rCellPopulation.GetNumNodes();
    mGradientsVegfAnalyticalApproxPde.resize(num_nodes, zero_vector<double>(DIM));

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin(); cell_iter != rCellPopulation.End(); ++cell_iter)
    {
        if (cell_iter->GetMutationState()->template IsType<VesselTipMutationState>())
        {
            // we collect the cell data necessary (node index and cell pointer)
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter); 

            c_vector<double, DIM> x_parent = rCellPopulation.GetLocationOfCellCentre((*cell_iter));
            c_vector<double, DIM> r_gradient_cell = zero_vector<double>(DIM);

            double Kc = sqrt((mDecayCoefficient-mCreationCoefficient)/mDiffusionCoefficient);
            double vegf_concentration_gradient = -Kc*(mSourceValue-mConstantBackground)*exp(-Kc*abs(x_parent[0]));
            
            r_gradient_cell(0) = mChiAnalyticalApproxPde*vegf_concentration_gradient; 

            mGradientsVegfAnalyticalApproxPde[node_index] = r_gradient_cell;
        }
    }
}

// Explicit instantiation
template class ChemoForceWithAnalyticalApproximationPde<1>;
template class ChemoForceWithAnalyticalApproximationPde<2>;
template class ChemoForceWithAnalyticalApproximationPde<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemoForceWithAnalyticalApproximationPde)