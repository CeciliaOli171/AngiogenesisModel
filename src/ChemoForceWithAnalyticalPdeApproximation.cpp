#include "ChemoForceWithAnalyticalPdeApproximation.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "VegfEquationPde.hpp"
#include "LinearBasisFunction.hpp"
#include "ReplicatableVector.hpp"

#include "Debug.hpp"

template<unsigned DIM>
ChemoForceWithAnalyticalPdeApproximation<DIM>::ChemoForceWithAnalyticalPdeApproximation(double chi,double diffusionCoefficient, double decayCoefficient, double creationCoefficient, double consumptionCoefficient, double boundaryCuboidMax, double maxValue)
    : ChemoForce<DIM>(chi, 0.0), mChiVegf(chi), mDiffusionCoefficient(diffusionCoefficient), mDecayCoefficient(decayCoefficient), mCreationCoefficient(creationCoefficient), mConsumptionCoefficient(consumptionCoefficient), mBoundaryCuboidMax(boundaryCuboidMax), mMaxValue(maxValue)
{
    assert(chi>0);
}

template<unsigned DIM>
ChemoForceWithAnalyticalPdeApproximation<DIM>::~ChemoForceWithAnalyticalPdeApproximation()
{
}

template<unsigned DIM>
c_vector<double, DIM>& ChemoForceWithAnalyticalPdeApproximation<DIM>::GetGradient(unsigned node_index)
{
    return mGradientsWithVegfApproximation[node_index];
}

template<unsigned DIM>
double ChemoForceWithAnalyticalPdeApproximation<DIM>::GetMagnitudeGradient(unsigned node_index)
{
    return norm_2(GetGradient(node_index))/mChiVegf;
}

template<unsigned DIM>
void ChemoForceWithAnalyticalPdeApproximation<DIM>::CalculateVegfGradient(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Initialise gradients size
    unsigned num_nodes = rCellPopulation.GetNumNodes();
    mGradientsWithVegfApproximation.resize(num_nodes, zero_vector<double>(DIM));

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin(); cell_iter != rCellPopulation.End(); ++cell_iter)
    {
        // we collect the cell data necessary (node index and cell pointer)
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        CellPtr pCell = rCellPopulation.GetCellUsingLocationIndex(node_index); 

        c_vector<double, DIM> x_parent = rCellPopulation.GetLocationOfCellCentre(pCell);

        c_vector<double, DIM> r_gradient_cell = zero_vector<double>(DIM);

        if (pCell->GetMutationState()->IsType<TipCellMutationState>())
        {
            double vegf_concentration_gradient = -sqrt((mDecayCoefficient-mCreationCoefficient)/mDiffusionCoefficient)*mMaxValue*(exp(-sqrt((mDecayCoefficient-mCreationCoefficient)/mDiffusionCoefficient)*x_parent[0]) + exp(sqrt((mDecayCoefficient-mCreationCoefficient)/mDiffusionCoefficient)*(x_parent[0]-2*mBoundaryCuboidMax)))/(1-exp(-sqrt(2*(mDecayCoefficient-mCreationCoefficient)/mDiffusionCoefficient)*mBoundaryCuboidMax));
            
            if(DIM == 3){
                r_gradient_cell(0) = mChiVegf*vegf_concentration_gradient; 
                r_gradient_cell(1) = 0.0; 
                r_gradient_cell(2) = 0.0; 
            } else if (DIM == 2){
                r_gradient_cell(0) = mChiVegf*vegf_concentration_gradient; 
                r_gradient_cell(1) = 0.0; 
            } else {
                r_gradient_cell(0) = mChiVegf*vegf_concentration_gradient; 
            }

            mGradientsWithVegfApproximation[node_index] += r_gradient_cell;
        }
    }
}

// Explicit instantiation
template class ChemoForceWithAnalyticalPdeApproximation<1>;
template class ChemoForceWithAnalyticalPdeApproximation<2>;
template class ChemoForceWithAnalyticalPdeApproximation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemoForceWithAnalyticalPdeApproximation)