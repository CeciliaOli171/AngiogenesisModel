#include "ConnectivityWriter.hpp"

#include "AbstractCellPopulation.hpp"
#include "CellLabel.hpp"
#include "CellData.hpp"

#include "Debug.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ConnectivityWriter<ELEMENT_DIM, SPACE_DIM>::ConnectivityWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("results.vizconnectivity")
{
    this->mVtkCellDataName = "Connectivity";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ConnectivityWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double connectivity = pCell->GetCellData()->GetItem("Connectivity");
    return connectivity;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConnectivityWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double connectivity = pCell->GetCellData()->GetItem("Connectivity");
    unsigned node_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
    *this->mpOutStream << node_index << " " << connectivity << " ";

    if(pCell->GetCellData()->GetItem("Anastomosis") != 0){
        unsigned connectivity_anastomosis = pCell->GetCellData()->GetItem("AnastomosisConnectivity");
        *this->mpOutStream << node_index << " " << connectivity_anastomosis << " ";
    }
}

// Explicit instantiation
template class ConnectivityWriter<1,1>;
template class ConnectivityWriter<1,2>;
template class ConnectivityWriter<2,2>;
template class ConnectivityWriter<1,3>;
template class ConnectivityWriter<2,3>;
template class ConnectivityWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ConnectivityWriter)