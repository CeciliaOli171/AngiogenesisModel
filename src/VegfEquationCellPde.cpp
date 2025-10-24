#include "VegfEquationCellPde.hpp"

#include "AbstractCellPopulation.hpp"
#include "TetrahedralMesh.hpp"
#include "StemCellProliferativeType.hpp"
#include "VesselTipMutationState.hpp"

#include "Debug.hpp"

template<unsigned DIM>
VegfEquationCellPde<DIM>::VegfEquationCellPde(AbstractCellPopulation<DIM,DIM>& rCellPopulation, double duDtCoefficient, double diffusionCoefficient, double decayCoefficient, double creationCoefficient, double consumptionCoefficient)
    : CellwiseSourceParabolicPde<DIM>(rCellPopulation, duDtCoefficient, diffusionCoefficient, decayCoefficient), mrCellPopulation(rCellPopulation), mDuDtCoefficient(duDtCoefficient), mDiffusionCoefficient(diffusionCoefficient), mDecayCoefficient(decayCoefficient), mCreationCoefficient(creationCoefficient), mConsumptionCoefficient(consumptionCoefficient)
{
}

template<unsigned DIM>
double VegfEquationCellPde<DIM>::GetDuDtCoefficient()
{
    return mDuDtCoefficient;
}

template<unsigned DIM>
double VegfEquationCellPde<DIM>::GetDiffusionCoefficient()
{
    return mDiffusionCoefficient;
}

template<unsigned DIM>
double VegfEquationCellPde<DIM>::GetDecayCoefficient()
{
    return mDecayCoefficient;
}

template<unsigned DIM>
double VegfEquationCellPde<DIM>::GetCreationCoefficient()
{
    return mCreationCoefficient;
}

template<unsigned DIM>
double VegfEquationCellPde<DIM>::GetConsumptionCoefficient()
{
    return mConsumptionCoefficient;
}

template<unsigned DIM>
double VegfEquationCellPde<DIM>::ComputeDuDtCoefficientFunction(const ChastePoint<DIM>& )
{
    return mDuDtCoefficient;
}

// LCOV_EXCL_START
template<unsigned DIM>
double VegfEquationCellPde<DIM>::ComputeSourceTerm(const ChastePoint<DIM>& rX, double u, Element<DIM,DIM>* pElement)
{
    NEVER_REACHED;
    return 0.0;
}
// LCOV_EXCL_STOP

template<unsigned DIM>
double VegfEquationCellPde<DIM>::ComputeSourceTermAtNode(const Node<DIM>& rNode, double u)
{
    double source_coefficient = 0.0;
    CellPtr pCell = mrCellPopulation.GetCellUsingLocationIndex(rNode.GetIndex());

    if (pCell->GetMutationState()->IsType<VesselTipMutationState>())
    {
        source_coefficient = mCreationCoefficient - mDecayCoefficient - mConsumptionCoefficient;
    }
    else 
    {
        source_coefficient = mCreationCoefficient - mDecayCoefficient;
    }

    // The source term is C*u
    return source_coefficient*u;
}

template<unsigned DIM>
c_matrix<double,DIM,DIM>  VegfEquationCellPde<DIM>::ComputeDiffusionTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement)
{
    return mDiffusionCoefficient*identity_matrix<double>(DIM);
}

// Explicit instantiation
template class VegfEquationCellPde<1>;
template class VegfEquationCellPde<2>;
template class VegfEquationCellPde<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VegfEquationCellPde)