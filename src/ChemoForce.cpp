#include "ChemoForce.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"

template<unsigned DIM>
ChemoForce<DIM>::ChemoForce(double chi, double cx, double cy)
    : AbstractForce<DIM>()
{
    assert(chi>0);
    assert(cx>0);
    mChi = chi;
    mCX = cx;
    mCY = cy;
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
double ChemoForce<DIM>::GetChemotacticGradientCoefficientYAxis()
{
    return mCY;
}

template<unsigned DIM>
void ChemoForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // initialisation 
    c_vector<double, DIM> chemoforce;
    c_vector<double,DIM> r_gradient;
    
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // we collect the cell data necessary (node index and cell pointer)
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        CellPtr pCell = rCellPopulation.GetCellUsingLocationIndex(node_index); 

        if (pCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
        {
            if(DIM == 3){
                r_gradient(0) = -mCX; 
                r_gradient(1) = -mCY; 
                r_gradient(2) = 0.0; 
            } else if (DIM == 2){
                r_gradient(0) = -mCX; 
                r_gradient(1) = -mCY; 
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
    *rParamsFile << "\t\t\t<Chemotactic Sensitivity>" << mChi << "</Chi>\n";
    *rParamsFile << "\t\t\t<Chemotactic Gradient Coefficient X Axis>" << mCX << "</CX>\n";
    *rParamsFile << "\t\t\t<Chemotactic Gradient Coefficient Y Axis>" << mCY << "</CY>\n";

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
