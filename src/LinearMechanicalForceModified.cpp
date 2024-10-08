#include "LinearMechanicalForceModified.hpp"
#include "CellLabel.hpp"
#include "CellData.hpp"

#include "BranchingCellMutationState.hpp"
#include "TipCellMutationState.hpp"
#include "VesselCellMutationState.hpp"

#include "Debug.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
LinearMechanicalForceModified<ELEMENT_DIM, SPACE_DIM>::LinearMechanicalForceModified()
   : GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
LinearMechanicalForceModified<ELEMENT_DIM, SPACE_DIM>::~LinearMechanicalForceModified()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double LinearMechanicalForceModified<ELEMENT_DIM, SPACE_DIM>::VariableSpringConstantMultiplicationFactor(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, bool isCloserThanRestLength)
{
    //TRACE("Begin Mechanical Force");

    // initialisation
    double ConstMultiFactor = 1.0;

    // we obtain the pointer of the cells considered 
    CellPtr pCellA = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
    CellPtr pCellB = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

    double branch_number_A = pCellA->GetCellData()->GetItem("BranchNumber");
    double branch_number_B = pCellB->GetCellData()->GetItem("BranchNumber");

    // if both the cells are on the same branch, then the factor is equal to 1
    // if the cells are not on the same branch, then they exert no attraction or repulsion to each other 
    if(branch_number_A == branch_number_B){
        ConstMultiFactor = 1.0;
    } else if (pCellA->GetMutationState() ->IsType<BranchingCellMutationState>()){
        ConstMultiFactor = 1.0;
    } else if (pCellB->GetMutationState() ->IsType<BranchingCellMutationState>()){
        ConstMultiFactor = 1.0;
    } else if (!isCloserThanRestLength)  {
        ConstMultiFactor = 0.0;
    }

    //TRACE("End Mechanical Force");

    return ConstMultiFactor;
}


// Explicit instantiation
template class LinearMechanicalForceModified<1,1>;
template class LinearMechanicalForceModified<1,2>;
template class LinearMechanicalForceModified<2,2>;
template class LinearMechanicalForceModified<1,3>;
template class LinearMechanicalForceModified<2,3>;
template class LinearMechanicalForceModified<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(LinearMechanicalForceModified)
