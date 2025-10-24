#include "ChemoForceWithAnalyticalApproximationPde.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "VegfEquationPde.hpp"
#include "LinearBasisFunction.hpp"
#include "ReplicatableVector.hpp"

#include "Debug.hpp"

template<unsigned DIM>
ChemoForceWithAnalyticalApproximationPde<DIM>::ChemoForceWithAnalyticalApproximationPde(double chiAnalyticalApproxPde, double hx, double diffusionCoefficient, double decayCoefficient, double creationCoefficient, double consumptionCoefficient, double sourceValue, double constantBackground)
    : ChemoForce<DIM>(chiAnalyticalApproxPde, hx, 0.0, 0.0, 0.0, 0.0, 0.0), mChiAnalyticalApproxPde(chiAnalyticalApproxPde), mHX(hx), mDiffusionCoefficient(diffusionCoefficient), mDecayCoefficient(decayCoefficient), mCreationCoefficient(creationCoefficient), mConsumptionCoefficient(consumptionCoefficient), mSourceValue(sourceValue), mConstantBackground(constantBackground)
{
    assert(chiAnalyticalApproxPde>0);
}

template<unsigned DIM>
ChemoForceWithAnalyticalApproximationPde<DIM>::~ChemoForceWithAnalyticalApproximationPde()
{
}

template<unsigned DIM>
c_vector<double, DIM>& ChemoForceWithAnalyticalApproximationPde<DIM>::GetGradient(AbstractCellPopulation<DIM>& rCellPopulation, unsigned node_index)
{
    CellPtr pCell = rCellPopulation.GetCellUsingLocationIndex(node_index);
    c_vector<double, DIM> x_parent = rCellPopulation.GetLocationOfCellCentre(pCell);
    c_vector<double, DIM> r_gradient_cell = zero_vector<double>(DIM);

    double Kc = sqrt((mDecayCoefficient-mCreationCoefficient)/mDiffusionCoefficient);
    double vegf_concentration_gradient = -Kc*(mSourceValue-mConstantBackground)*exp(-Kc*abs(x_parent[0]));
    
    r_gradient_cell(0) = mChiAnalyticalApproxPde*vegf_concentration_gradient-mHX; 

    return r_gradient_cell;
}

// Explicit instantiation
template class ChemoForceWithAnalyticalApproximationPde<1>;
template class ChemoForceWithAnalyticalApproximationPde<2>;
template class ChemoForceWithAnalyticalApproximationPde<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemoForceWithAnalyticalApproximationPde)