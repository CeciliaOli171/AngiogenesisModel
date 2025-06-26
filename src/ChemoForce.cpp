#include "ChemoForce.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "VegfEquationPde.hpp"
#include "LinearBasisFunction.hpp"
#include "ReplicatableVector.hpp"

#include "Debug.hpp"

template<unsigned DIM>
ChemoForce<DIM>::ChemoForce(double chi, double cx)
    : AbstractForce<DIM>(), mChi(chi)
{
    assert(chi>0);
    assert(cx>0);
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
double ChemoForce<DIM>::GetChemotacticGradientCoefficientXAxis()
{
    return mCX;
}

template<unsigned DIM>
c_vector<double, DIM>& ChemoForce<DIM>::GetGradient(unsigned node_index)
{
    return mGradients[node_index];
}

template<unsigned DIM>
double ChemoForce<DIM>::GetMagnitudeGradient(unsigned node_index)
{
    return norm_2(GetGradient(node_index))/mChi;
}

template<unsigned DIM>
void ChemoForce<DIM>::CalculateVegfGradient(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Initialise gradients size
    unsigned num_nodes = rCellPopulation.GetNumNodes();
    mGradients.resize(num_nodes, zero_vector<double>(DIM));

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin(); cell_iter != rCellPopulation.End(); ++cell_iter)
    {
        // we collect the cell data necessary (node index and cell pointer)
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        CellPtr pCell = rCellPopulation.GetCellUsingLocationIndex(node_index); 

        c_vector<double, DIM> r_gradient_cell = zero_vector<double>(DIM);

        if (pCell->GetMutationState()->IsType<TipCellMutationState>())
        {
            if(DIM == 3){
                r_gradient_cell(0) = -*mCX; 
                r_gradient_cell(1) = 0.0; 
                r_gradient_cell(2) = 0.0; 
            } else if (DIM == 2){
                r_gradient_cell(0) = -*mCX; 
                r_gradient_cell(1) = 0.0; 
            } else {
                r_gradient_cell(0) = -*mCX; 
            }

            mGradients[node_index] += r_gradient_cell;
        }
    }
}

template<unsigned DIM>
void ChemoForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    //TRACE("Begin Chemotactic Force");

    // initialisation 
    c_vector<double, DIM> chemoforce;
    c_vector<double, DIM> r_gradient;

    // we calculate the gradient of the solution of the vegf pde 
    CalculateVegfGradient(rCellPopulation);
    
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // we collect the cell data necessary (node index and cell pointer)
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        CellPtr pCell = rCellPopulation.GetCellUsingLocationIndex(node_index); 

        if (pCell->GetMutationState()->IsType<TipCellMutationState>())
        {
            // we collect the gradient at the cell position 
            c_vector<double, DIM> r_gradient_cell = GetGradient(node_index);
            double magnitude_gradient_cell = GetMagnitudeGradient(node_index);

            // force += chi * gradC
            if(magnitude_gradient_cell != 0.0){
                // chemoforce = r_gradient_cell/magnitude_gradient_cell;
                chemoforce = mChi*r_gradient_cell;
            } else {
                chemoforce = zero_vector<double>(DIM);
            }
            rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(chemoforce);
        }
    }

    //TRACE("End Chemotactic Force");
}

template<unsigned DIM>
void ChemoForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<Chemotactic Sensitivity>" << mChi << "</Chi>\n";

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
