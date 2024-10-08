#include "BranchesNumberWriter.hpp"

#include "AbstractCellPopulation.hpp"
#include "CellLabel.hpp"
#include "CellData.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
BranchesNumberWriter<ELEMENT_DIM, SPACE_DIM>::BranchesNumberWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("results.vizbranchesnumber")
{
    this->mVtkCellDataName = "Branches Number";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double BranchesNumberWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double branches_number = pCell->GetCellData()->GetItem("BranchNumber");
    return branches_number;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BranchesNumberWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double branches_number = pCell->GetCellData()->GetItem("BranchNumber");
    *this->mpOutStream << branches_number << " ";
}

// Explicit instantiation
template class BranchesNumberWriter<1,1>;
template class BranchesNumberWriter<1,2>;
template class BranchesNumberWriter<2,2>;
template class BranchesNumberWriter<1,3>;
template class BranchesNumberWriter<2,3>;
template class BranchesNumberWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(BranchesNumberWriter)