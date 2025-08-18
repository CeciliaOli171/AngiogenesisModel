#include "PersistenceForce.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "CellData.hpp"

#include "SimulationTime.hpp"
#include "Debug.hpp"

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
    //TRACE("Begin Persistence Force");
    
    for(typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter) 
    {        
        if (cell_iter->GetMutationState()->template IsType<TipCellMutationState>()) // && dt % 180 == 0 
        {
            c_vector<double, DIM> persistenceforce = zero_vector<double>(DIM); 

            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

            std::vector<std::string> mDataLabels = std::vector<std::string>{ "fx", "fy", "fz" };
            for(unsigned i = 0; i < DIM; ++i){
                double f = (*cell_iter)->GetCellData()->GetItem(mDataLabels[i]);
                //PRINT_VARIABLE(f)
                persistenceforce(i) = mOmegap*f;
            
            }
            
            rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(persistenceforce);

            // test time force 
            //for(int i = 0; i < 100; ++i){rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(persistenceforce);}
        }

    }

    //TRACE("End Persistence Force");
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