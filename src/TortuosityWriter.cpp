#include "TortuosityWriter.hpp"

#include "AbstractCellPopulation.hpp"
#include "CellLabel.hpp"
#include "CellData.hpp"

#include "Debug.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
TortuosityWriter<ELEMENT_DIM, SPACE_DIM>::TortuosityWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("results.viztortuosity")
{
    this->mVtkCellDataName = "Tortuosity";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double TortuosityWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double tortuosity = pCell->GetCellData()->GetItem("TortuosityParent");
    return tortuosity;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TortuosityWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double tortuosity = pCell->GetCellData()->GetItem("TortuosityParent");
    const c_vector<double, SPACE_DIM> tipcell_coordinates = pCellPopulation->GetLocationOfCellCentre(pCell);

    if(SPACE_DIM == 1){
        *this->mpOutStream << tortuosity << " ";
        *this->mpOutStream << tipcell_coordinates(0) << " ";
    } else if(SPACE_DIM == 2){
        *this->mpOutStream << tortuosity << " ";
        *this->mpOutStream << tipcell_coordinates(0) << " ";
        *this->mpOutStream << tipcell_coordinates(1) << " ";
    } else if(SPACE_DIM == 3){
        *this->mpOutStream << tortuosity << " ";
        *this->mpOutStream << tipcell_coordinates(0) << " ";
        *this->mpOutStream << tipcell_coordinates(1) << " ";
        *this->mpOutStream << tipcell_coordinates(2) << " ";
    }

}

// Explicit instantiation
template class TortuosityWriter<1,1>;
template class TortuosityWriter<1,2>;
template class TortuosityWriter<2,2>;
template class TortuosityWriter<1,3>;
template class TortuosityWriter<2,3>;
template class TortuosityWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(TortuosityWriter)