#include "ChemoForceTwoLesionsYAxis.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "VegfEquationPde.hpp"
#include "LinearBasisFunction.hpp"
#include "ReplicatableVector.hpp"

#include "Debug.hpp"

template<unsigned DIM>
ChemoForceTwoLesionsYAxis<DIM>::ChemoForceTwoLesionsYAxis(double chiTwoLesionsYAxis, double hy, double diffusionCoefficient, double decayCoefficient, double creationCoefficient, double consumptionCoefficient, double sourceValue1, double sourceValue2, double constantBackground, double locationLesion1, double locationLesion2)
    : ChemoForce<DIM>(chiTwoLesionsYAxis, 0.0, hy, 0.0, 0.0, 0.0, 0.0), mChiTwoLesionsYAxis(chiTwoLesionsYAxis), mDiffusionCoefficient(diffusionCoefficient), mDecayCoefficient(decayCoefficient), mCreationCoefficient(creationCoefficient), mConsumptionCoefficient(consumptionCoefficient), mSourceValue1(sourceValue1), mSourceValue2(sourceValue2), mConstantBackground(constantBackground), mLocationLesion1(locationLesion1), mLocationLesion2(locationLesion2)
{
    assert(chiTwoLesionsYAxis>0);
}

template<unsigned DIM>
ChemoForceTwoLesionsYAxis<DIM>::~ChemoForceTwoLesionsYAxis()
{
}

template<unsigned DIM>
c_vector<double, DIM>& ChemoForceTwoLesionsYAxis<DIM>::GetGradient(unsigned node_index)
{
    return mGradientsVegfTwoLesionsYAxis[node_index];
}

template<unsigned DIM>
double ChemoForceTwoLesionsYAxis<DIM>::GetMagnitudeGradient(unsigned node_index)
{
    return norm_2(GetGradient(node_index));
}

template<unsigned DIM>
void ChemoForceTwoLesionsYAxis<DIM>::CalculateVegfGradient(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Initialise gradients size
    unsigned num_nodes = rCellPopulation.GetNumNodes();
    mGradientsVegfTwoLesionsYAxis.resize(num_nodes, zero_vector<double>(DIM));

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin(); cell_iter != rCellPopulation.End(); ++cell_iter)
    {
        if (cell_iter->GetMutationState()->template IsType<VesselTipMutationState>())
        {
            // we collect the cell data necessary (node index and cell pointer)
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

            c_vector<double, DIM> x_parent = rCellPopulation.GetLocationOfCellCentre((*cell_iter));
            c_vector<double, DIM> r_gradient_cell = zero_vector<double>(DIM);

            double Kc = sqrt((mDecayCoefficient-mCreationCoefficient)/mDiffusionCoefficient);
            double vegf_concentration_gradient_1 = -Kc*(mSourceValue1-mConstantBackground)*exp(-Kc*(x_parent[0]-mLocationLesion1));
            double vegf_concentration_gradient_2 = Kc*(mSourceValue2-mConstantBackground)*exp(Kc*(x_parent[0]-mLocationLesion2));

            double vegf_concentration_gradient = (vegf_concentration_gradient_1 + vegf_concentration_gradient_2)/2;
            
            r_gradient_cell(0) = mChiTwoLesionsYAxis*vegf_concentration_gradient; 

            mGradientsVegfTwoLesionsYAxis[node_index] = r_gradient_cell;
        }
    }
}

// Explicit instantiation
template class ChemoForceTwoLesionsYAxis<1>;
template class ChemoForceTwoLesionsYAxis<2>;
template class ChemoForceTwoLesionsYAxis<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemoForceTwoLesionsYAxis)