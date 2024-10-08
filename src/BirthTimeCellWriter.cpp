#include "BirthTimeCellWriter.hpp"

#include "AbstractCellPopulation.hpp"
#include "CellLabel.hpp"
#include "CellData.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
BirthTimeCellWriter<ELEMENT_DIM, SPACE_DIM>::BirthTimeCellWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("results.vizbirthtimecell")
{
    this->mVtkCellDataName = "Birth Time Cell";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double BirthTimeCellWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double birth_time_cell = pCell->GetCellData()->GetItem("BirthTimeCell");
    return birth_time_cell;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BirthTimeCellWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double birth_time_cell = pCell->GetCellData()->GetItem("BirthTimeCell");
    *this->mpOutStream << birth_time_cell << " ";
}

// Explicit instantiation
template class BirthTimeCellWriter<1,1>;
template class BirthTimeCellWriter<1,2>;
template class BirthTimeCellWriter<2,2>;
template class BirthTimeCellWriter<1,3>;
template class BirthTimeCellWriter<2,3>;
template class BirthTimeCellWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(BirthTimeCellWriter)