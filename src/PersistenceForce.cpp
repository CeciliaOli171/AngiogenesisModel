#include "PersistenceForce.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "CellData.hpp"

#include "SimulationTime.hpp"

template<unsigned DIM>
PersistenceForce<DIM>::PersistenceForce(double omegap)
    : AbstractForce<DIM>()
{
    assert(omegap > 0);
    mOmegap = omegap;
}

template<unsigned DIM>
PersistenceForce<DIM>::~PersistenceForce()
{
}

template<unsigned DIM>
double PersistenceForce<DIM>::GetPersistenceCoefficient()
{
    return mOmegap;
}

// overrides AddForceContribution()
template<unsigned DIM>
void PersistenceForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // initialisation 
    c_vector<double, DIM> persistenceforce; 

    int dt = SimulationTime::Instance()->GetTime();
    
    for(typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter) 
    {
        // we collect the cell data (node index and cell pointer)
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        CellPtr pCell = rCellPopulation.GetCellUsingLocationIndex(node_index);

        if (pCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()) // && dt % 180 == 0 
        {
            // we find back the location of the cell at the previous time step by using CellData
            c_vector<double,DIM> old_r_cellmovement;
            double old_x_coordinate = cell_iter->GetCellData()->GetItem("old_x_coordinate");  
            old_r_cellmovement(0) = old_x_coordinate; 
            if(DIM > 1){
                double old_y_coordinate = cell_iter->GetCellData()->GetItem("old_y_coordinate");  
                old_r_cellmovement(1) = old_y_coordinate;
            } 
            if(DIM == 3){
                double old_z_coordinate = cell_iter->GetCellData()->GetItem("old_z_coordinate");  
                old_r_cellmovement(2) = old_z_coordinate;
            }
        
            c_vector<double, DIM> new_r_cellmovement = rCellPopulation.GetLocationOfCellCentre(*cell_iter); 
            c_vector<double,DIM> r_cellmovement = new_r_cellmovement - old_r_cellmovement; 
            double magnitude_of_cellmovement = norm_2(r_cellmovement);

            if(magnitude_of_cellmovement != 0.0){
                persistenceforce = dt*(mOmegap/magnitude_of_cellmovement)*r_cellmovement;
            } else {
                persistenceforce = zero_vector<double>(DIM);
            }
            
            rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(persistenceforce);
            
            // we store in cell data the position of the centre of the cells to use them later 
            // PBLM: only store the previous persistence force coordinates and not the ones with all the forces applied 
            // TO CHANGE 
            cell_iter->GetCellData()->SetItem("old_x_coordinate", new_r_cellmovement(0)); 

            if(DIM < 1){
                cell_iter->GetCellData()->SetItem("old_y_coordinate", new_r_cellmovement(1));  
            } 
            if(DIM == 3){
                cell_iter->GetCellData()->SetItem("old_z_coordinate", new_r_cellmovement(2));  
            }
        }

    }
}

template<unsigned DIM>
void PersistenceForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<Directional Persistence Coefficient>" << mOmegap << "</Omegap>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class PersistenceForce<1>;
template class PersistenceForce<2>;
template class PersistenceForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PersistenceForce)