#include "AnastomosisWriter.hpp"

#include "AbstractCellPopulation.hpp"
#include "CellLabel.hpp"
#include "CellData.hpp"

#include "Debug.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AnastomosisWriter<ELEMENT_DIM, SPACE_DIM>::AnastomosisWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("results.vizanastomosis")
{
    this->mVtkCellDataName = "anastomosis";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AnastomosisWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double anastomosis = pCell->GetCellData()->GetItem("Anastomosis");
    return anastomosis;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AnastomosisWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double anastomosis = pCell->GetCellData()->GetItem("Anastomosis");
    *this->mpOutStream << anastomosis << " ";

}

// Explicit instantiation
template class AnastomosisWriter<1,1>;
template class AnastomosisWriter<1,2>;
template class AnastomosisWriter<2,2>;
template class AnastomosisWriter<1,3>;
template class AnastomosisWriter<2,3>;
template class AnastomosisWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(AnastomosisWriter)