#include "ChemoForce.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"

template<unsigned DIM>
ChemoForce<DIM>::ChemoForce(double chi, double cx)
    : AbstractForce<DIM>()
{
    assert(chi>0);
    assert(cx>0);
    mChi = chi;
    mCX = cx;
}

template<unsigned DIM>
ChemoForce<DIM>::~ChemoForce()
{
}

template<unsigned DIM>
double ChemoForce<DIM>::GetChemotacticSensitivity()
{
    return mChi;
}

template<unsigned DIM>
void ChemoForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // we create the vector force 
    c_vector<double, DIM> chemoforce = zero_vector<double>(DIM);
    
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        CellPtr cell_ptr = rCellPopulation.GetCellUsingLocationIndex(node_index); 
        if (cell_ptr->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
        {
            c_vector<double,DIM> r_gradient;
            if(DIM == 3){
                r_gradient(0) = -mCX; 
                r_gradient(1) = 0.0; 
                r_gradient(2) = 0.0; 
            } else if (DIM == 2){
                r_gradient(0) = -mCX; 
                r_gradient(1) = 0.0; 
            } else {
                r_gradient(0) = -mCX; 
            }

            // force += chi * gradC
            chemoforce = mChi*r_gradient;
            rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(chemoforce);
        }
    }
}

template<unsigned DIM>
void ChemoForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // No parameters to include

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class ChemoForce<1>;
template class ChemoForce<2>;
template class ChemoForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemoForce)
