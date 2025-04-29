#include "ChemoForceWithAnalyticalApproximationPde.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "VegfEquationPde.hpp"
#include "LinearBasisFunction.hpp"
#include "ReplicatableVector.hpp"

#include "Debug.hpp"

template<unsigned DIM>
ChemoForceWithAnalyticalApproximationPde<DIM>::ChemoForceWithAnalyticalApproximationPde(double chiAnalyticalApproxPde, double diffusionCoefficient, double decayCoefficient, double creationCoefficient, double consumptionCoefficient, double sourceValue, double constantBackground, double boundaryCuboidMax, double maxValue)
    : ChemoForce<DIM>(chiAnalyticalApproxPde, 0.0), mChiAnalyticalApproxPde(chiAnalyticalApproxPde), mDiffusionCoefficient(diffusionCoefficient), mDecayCoefficient(decayCoefficient), mCreationCoefficient(creationCoefficient), mConsumptionCoefficient(consumptionCoefficient), mSourceValue(sourceValue), mConstantBackground(constantBackground), mBoundaryCuboidMax(boundaryCuboidMax), mMaxValue(maxValue)
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
    return norm_2(GetGradient(node_index))/mChiAnalyticalApproxPde;
}

template<unsigned DIM>
void ChemoForceWithAnalyticalApproximationPde<DIM>::CalculateVegfGradient(AbstractCellPopulation<DIM>& rCellPopulation)
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
            // double vegf_concentration_gradient = -sqrt((mDecayCoefficient-mCreationCoefficient)/mDiffusionCoefficient)*mMaxValue*(exp(-sqrt((mDecayCoefficient-mCreationCoefficient)/mDiffusionCoefficient)*x_parent[0]) + exp(sqrt((mDecayCoefficient-mCreationCoefficient)/mDiffusionCoefficient)*(x_parent[0]-2*mBoundaryCuboidMax)))/(1-exp(-sqrt(2*(mDecayCoefficient-mCreationCoefficient)/mDiffusionCoefficient)*mBoundaryCuboidMax)); // old version

            double Kc = sqrt((mDecayCoefficient-mCreationCoefficient)/mDiffusionCoefficient);
            double vegf_concentration_gradient = -Kc*(mSourceValue-mConstantBackground)*exp(-Kc*x_parent[0]);
            
            if(DIM == 3){
                r_gradient_cell(0) = mChiAnalyticalApproxPde*vegf_concentration_gradient; 
                r_gradient_cell(1) = 0.0; 
                r_gradient_cell(2) = 0.0; 
            } else if (DIM == 2){
                r_gradient_cell(0) = mChiAnalyticalApproxPde*vegf_concentration_gradient; 
                r_gradient_cell(1) = 0.0; 
            } else {
                r_gradient_cell(0) = mChiAnalyticalApproxPde*vegf_concentration_gradient; 
            }

            mGradientsVegfAnalyticalApproxPde[node_index] += r_gradient_cell;
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