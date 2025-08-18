#include "ChemoForceTwoLesionsXAxis.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "VegfEquationPde.hpp"
#include "LinearBasisFunction.hpp"
#include "ReplicatableVector.hpp"

#include "Debug.hpp"

template<unsigned DIM>
ChemoForceTwoLesionsXAxis<DIM>::ChemoForceTwoLesionsXAxis(double chiTwoLesionsXAxis, double hx, double diffusionCoefficient, double decayCoefficient, double creationCoefficient, double consumptionCoefficient, double sourceValue1, double sourceValue2, double constantBackground, double locationLesion1, double locationLesion2)
    : ChemoForce<DIM>(chiTwoLesionsXAxis, 1E-4, hx, 0.0, 0.0, 0.0, 0.0), mChiTwoLesionsXAxis(chiTwoLesionsXAxis), mDiffusionCoefficient(diffusionCoefficient), mDecayCoefficient(decayCoefficient), mCreationCoefficient(creationCoefficient), mConsumptionCoefficient(consumptionCoefficient), mSourceValue1(sourceValue1), mSourceValue2(sourceValue2), mConstantBackground(constantBackground), mLocationLesion1(locationLesion1), mLocationLesion2(locationLesion2)
{
    assert(chiTwoLesionsXAxis>0);
}

template<unsigned DIM>
ChemoForceTwoLesionsXAxis<DIM>::~ChemoForceTwoLesionsXAxis()
{
}

template<unsigned DIM>
c_vector<double, DIM>& ChemoForceTwoLesionsXAxis<DIM>::GetGradient(unsigned node_index)
{
    return mGradientsVegfTwoLesionsXAxis[node_index];
}

template<unsigned DIM>
double ChemoForceTwoLesionsXAxis<DIM>::GetMagnitudeGradient(unsigned node_index)
{
    return norm_2(GetGradient(node_index));
}

template<unsigned DIM>
void ChemoForceTwoLesionsXAxis<DIM>::CalculateVegfGradient(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Initialise gradients size
    unsigned num_nodes = rCellPopulation.GetNumNodes();
    mGradientsVegfTwoLesionsXAxis.resize(num_nodes, zero_vector<double>(DIM));

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin(); cell_iter != rCellPopulation.End(); ++cell_iter)
    {
        if (cell_iter->GetMutationState()->template IsType<TipCellMutationState>())
        {
            // we collect the cell data necessary (node index and cell pointer)
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

            c_vector<double, DIM> x_parent = rCellPopulation.GetLocationOfCellCentre((*cell_iter));
            c_vector<double, DIM> r_gradient_cell = zero_vector<double>(DIM);

            double Kc = sqrt((mDecayCoefficient-mCreationCoefficient)/mDiffusionCoefficient);
            double vegf_concentration_gradient_1 = -Kc*(mSourceValue1-mConstantBackground)*exp(-Kc*(x_parent[0]-mLocationLesion1));
            double vegf_concentration_gradient_2 = -Kc*(mSourceValue2-mConstantBackground)*exp(-Kc*(x_parent[0]-mLocationLesion2));
            double vegf_concentration_gradient = (vegf_concentration_gradient_1+vegf_concentration_gradient_2)/2;
            
            r_gradient_cell(0) = mChiTwoLesionsXAxis*vegf_concentration_gradient; 
            
            mGradientsVegfTwoLesionsXAxis[node_index] = r_gradient_cell;
        }
    }
}

// Explicit instantiation
template class ChemoForceTwoLesionsXAxis<1>;
template class ChemoForceTwoLesionsXAxis<2>;
template class ChemoForceTwoLesionsXAxis<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemoForceTwoLesionsXAxis)