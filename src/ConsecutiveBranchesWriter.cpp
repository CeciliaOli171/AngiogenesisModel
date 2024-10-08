#include "ConsecutiveBranchesWriter.hpp"

#include "AbstractCellPopulation.hpp"
#include "CellLabel.hpp"
#include "CellData.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ConsecutiveBranchesWriter<ELEMENT_DIM, SPACE_DIM>::ConsecutiveBranchesWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("results.vizconsecutivebranches")
{
    this->mVtkCellDataName = "Consecutive Branches";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ConsecutiveBranchesWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double consecutive_branch_number = pCell->GetCellData()->GetItem("ConsecutiveBranches");
    return consecutive_branch_number;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConsecutiveBranchesWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double consecutive_branch_number = pCell->GetCellData()->GetItem("ConsecutiveBranches");
    *this->mpOutStream << consecutive_branch_number << " ";
}

// Explicit instantiation
template class ConsecutiveBranchesWriter<1,1>;
template class ConsecutiveBranchesWriter<1,2>;
template class ConsecutiveBranchesWriter<2,2>;
template class ConsecutiveBranchesWriter<1,3>;
template class ConsecutiveBranchesWriter<2,3>;
template class ConsecutiveBranchesWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ConsecutiveBranchesWriter)