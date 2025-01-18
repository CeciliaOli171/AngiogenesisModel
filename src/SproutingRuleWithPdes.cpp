#include "SproutingRuleWithPdes.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "VesselCellMutationState.hpp"
#include "TipCellMutationState.hpp"
#include "BranchingCellMutationState.hpp"
#include "UniformCellCycleModel.hpp"

#include "RandomNumberGenerator.hpp"
#include <algorithm>
#include "Debug.hpp"
#include <random>
#include <cstdlib>
#include <iterator>
#include "ReplicatableVector.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SproutingRuleWithPdes<ELEMENT_DIM, SPACE_DIM>::SproutingRuleWithPdes(double MaxSproutingRate, boost::shared_ptr<AbstractBoxDomainPdeModifier<SPACE_DIM> > pPdeModifier)
    : SproutingRule<ELEMENT_DIM, SPACE_DIM>(MaxSproutingRate), mpPdeModifier(pPdeModifier)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SproutingRuleWithPdes<ELEMENT_DIM, SPACE_DIM>::~SproutingRuleWithPdes()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double SproutingRuleWithPdes<ELEMENT_DIM, SPACE_DIM>::GetSproutingProbability(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pParentCell){
    // initialisation 
    double Psprout;

    // we need to access the vegf concentration at the cell position x 
    TetrahedralMesh<SPACE_DIM,SPACE_DIM>* p_coarse_mesh = mpPdeModifier->GetFeMesh();
    ReplicatableVector previous_solution(mpPdeModifier->GetSolution());
    unsigned elem_index = mpPdeModifier->GetFeMesh()->GetContainingElementIndex(rCellPopulation.GetLocationOfCellCentre(pParentCell));
    Element<SPACE_DIM,SPACE_DIM>* p_element = p_coarse_mesh->GetElement(elem_index);
    double vegf_concentration = previous_solution[p_element->GetNodeGlobalIndex(0)];

    Psprout = mMaxSproutingRate*pow(vegf_concentration,1); // test since the concentration is between 0 and 1 

    return Psprout;
}

// Explicit instantiation
template class SproutingRuleWithPdes<1,1>;
template class SproutingRuleWithPdes<1,2>;
template class SproutingRuleWithPdes<1,3>;
template class SproutingRuleWithPdes<2,2>;
template class SproutingRuleWithPdes<2,3>;
template class SproutingRuleWithPdes<3,3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SproutingRuleWithPdes)