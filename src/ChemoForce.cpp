#include "ChemoForce.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "VegfEquationPde.hpp"
#include "LinearBasisFunction.hpp"
#include "ReplicatableVector.hpp"

#include "Debug.hpp"

template<unsigned DIM>
ChemoForce<DIM>::ChemoForce(double chi, double hx, double hy, double hz, double cx, double cy, double cz)
    : AbstractForce<DIM>(), mChi(chi), mHX(hx), mHY(hy), mHZ(hz), mCX(cx), mCY(cy), mCZ(cz)
{
    assert(chi>0);

    assert(hx>0);
    assert(hy>0);
    assert(hz>0);
    mHX = hx;
    mHY = hy;
    mHZ = hz;

    assert(cx>0);
    assert(cy>0);
    assert(cz>0);
    mCX = cx;
    mCY = cy;
    mCZ = cz;
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
double ChemoForce<DIM>::GetChemoattractantGradientFactorXAxis()
{
    return mHX;
}

template<unsigned DIM>
double ChemoForce<DIM>::GetChemoattractantGradientFactorYAxis()
{
    return mHY;
}

template<unsigned DIM>
double ChemoForce<DIM>::GetChemoattractantGradientFactorZAxis()
{
    return mHZ;
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
double ChemoForce<DIM>::GetChemotacticGradientCoefficientZAxis()
{
    return mCZ;
}

template<unsigned DIM>
c_vector<double, DIM>& ChemoForce<DIM>::GetGradient(unsigned node_index)
{
    return mGradients[node_index];
}

template<unsigned DIM>
double ChemoForce<DIM>::GetMagnitudeGradient(unsigned node_index)
{
    return norm_2(GetGradient(node_index));
}

template<unsigned DIM>
void ChemoForce<DIM>::CalculateVegfGradient(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Initialise gradients size
    unsigned num_nodes = rCellPopulation.GetNumNodes();
    mGradients.resize(num_nodes, zero_vector<double>(DIM));

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin(); cell_iter != rCellPopulation.End(); ++cell_iter)
    {
        if ((*cell_iter)->GetMutationState()->template IsType<VesselTipMutationState>())
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

            c_vector<double, DIM> r_gradient_cell = zero_vector<double>(DIM);

            if(DIM == 3){
                r_gradient_cell(0) = -mCX; 
                r_gradient_cell(1) = -mCY; 
                r_gradient_cell(2) = -mCZ; 
            } else if (DIM == 2){
                r_gradient_cell(0) = -mCX; 
                r_gradient_cell(1) = -mCY; 
            } else {
                r_gradient_cell(0) = -mCX; 
            }

            mGradients[node_index] = r_gradient_cell;
        }
    }
}

template<unsigned DIM>
void ChemoForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // we calculate the gradient of the solution of the vegf pde 
    CalculateVegfGradient(rCellPopulation);
    
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        if ((*cell_iter)->GetMutationState()->template IsType<VesselTipMutationState>()) 
        {
            // initialisation 
            c_vector<double, DIM> chemoforce = zero_vector<double>(DIM);

            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

            // we collect the gradient at the cell position 
            c_vector<double, DIM> r_gradient_cell = GetGradient(node_index);
            double magnitude_gradient_cell = GetMagnitudeGradient(node_index);

            if(DIM == 3){
                r_gradient_cell(0) = r_gradient_cell(0)-mHX; 
                r_gradient_cell(1) = r_gradient_cell(1)-mHY; 
                r_gradient_cell(2) = r_gradient_cell(2)-mHZ; 
            } else if (DIM == 2){
                r_gradient_cell(0) = r_gradient_cell(0)-mHX; 
                r_gradient_cell(1) = r_gradient_cell(1)-mHY; 
            } else {
                r_gradient_cell(0) = r_gradient_cell(0)-mHX; 
            }

            // force += chi * gradC
            if(magnitude_gradient_cell != 0.0){
                chemoforce = mChi*r_gradient_cell;
            } 

            rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(chemoforce);

            // test time force
            // for(int i = 0; i < 100; ++i){rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(chemoforce);}
        }
    }
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
