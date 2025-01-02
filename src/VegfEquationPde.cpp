#include "VegfEquationPde.hpp"

#include "ApoptoticCellProperty.hpp"
#include "AbstractCellPopulation.hpp"
#include "TetrahedralMesh.hpp"
#include "StemCellProliferativeType.hpp"
#include "TipCellMutationState.hpp"

#include "Debug.hpp"

template<unsigned DIM>
VegfEquationPde<DIM>::VegfEquationPde(AbstractCellPopulation<DIM,DIM>& rCellPopulation,
                                                            double duDtCoefficient,
                                                            double diffusionCoefficient,
                                                            double sourceCoefficient,
                                                            double creationCoefficient,
                                                            double consumptionCoefficient)
    : AveragedSourceParabolicPde<DIM>(rCellPopulation, duDtCoefficient, diffusionCoefficient, sourceCoefficient),
      mrCellPopulation(rCellPopulation),
      mDuDtCoefficient(duDtCoefficient),
      mDiffusionCoefficient(diffusionCoefficient),
      mSourceCoefficient(sourceCoefficient),
      mCreationCoefficient(creationCoefficient),
      mConsumptionCoefficient(consumptionCoefficient)
{
}

// template<unsigned DIM>
// VegfEquationPde<DIM>:: ~VegfEquationPde(){
// }

template<unsigned DIM>
void VegfEquationPde<DIM>::SetupSourceTerms(TetrahedralMesh<DIM,DIM>& rCoarseMesh, std::map<CellPtr, unsigned>* pCellPdeElementMap) // must be called before solve
{
    // Allocate memory
    mCellDensityOnCoarseElements.resize(rCoarseMesh.GetNumElements());
    for (unsigned elem_index=0; elem_index<mCellDensityOnCoarseElements.size(); elem_index++)
    {
        mCellDensityOnCoarseElements[elem_index] = 0.0;
    }

    // Loop over cells, find which coarse element it is in, and add 1 to mSourceTermOnCoarseElements[elem_index]
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = mrCellPopulation.Begin();
         cell_iter != mrCellPopulation.End();
         ++cell_iter)
    {
        // test for angiogenesis model pde 
        // we only consider tip cells for the vegf concentration PDE 
        // we collect the cell data necessary (node index and cell pointer)
        unsigned node_index = mrCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        CellPtr pCell = mrCellPopulation.GetCellUsingLocationIndex(node_index);

        unsigned elem_index = 0;
        const ChastePoint<DIM>& r_position_of_cell = mrCellPopulation.GetLocationOfCellCentre(*cell_iter);

        if (pCellPdeElementMap != nullptr)
        {
            elem_index = (*pCellPdeElementMap)[*cell_iter];
        }
        else
        {
            elem_index = rCoarseMesh.GetContainingElementIndex(r_position_of_cell);
        }

        //test for angiogenesis model pde
        if (pCell->GetCellProliferativeType()->IsType<TipCellMutationState>()){
            mCellDensityOnCoarseElements[elem_index] += 1.0;
        }
    }

    // Then divide each entry of mSourceTermOnCoarseElements by the element's area
    c_matrix<double, DIM, DIM> jacobian;
    double det;
    for (unsigned elem_index=0; elem_index<mCellDensityOnCoarseElements.size(); elem_index++)
    {
        rCoarseMesh.GetElement(elem_index)->CalculateJacobian(jacobian, det);
        mCellDensityOnCoarseElements[elem_index] /= rCoarseMesh.GetElement(elem_index)->GetVolume(det);
    }
}

template<unsigned DIM>
double VegfEquationPde<DIM>::ComputeDuDtCoefficientFunction(const ChastePoint<DIM>& )
{
    return mDuDtCoefficient;
}

template<unsigned DIM>
double VegfEquationPde<DIM>::ComputeSourceTerm(const ChastePoint<DIM>& rX, double u, Element<DIM,DIM>* pElement)
{
    assert(!mCellDensityOnCoarseElements.empty());
    double coefficient = mCreationCoefficient - mConsumptionCoefficient * mCellDensityOnCoarseElements[pElement->GetIndex()];

    // The source term is C*u
    return coefficient*u;
}

// LCOV_EXCL_START
template<unsigned DIM>
double VegfEquationPde<DIM>::ComputeSourceTermAtNode(const Node<DIM>& rNode, double u)
{
    NEVER_REACHED;
    return 0.0;
}
// LCOV_EXCL_STOP

template<unsigned DIM>
c_matrix<double,DIM,DIM> VegfEquationPde<DIM>::ComputeDiffusionTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement)
{
    return mDiffusionCoefficient*identity_matrix<double>(DIM);
}

// Explicit instantiation
template class VegfEquationPde<1>;
template class VegfEquationPde<2>;
template class VegfEquationPde<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VegfEquationPde)
