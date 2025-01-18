#include "SproutingRuleWithAnalyticalPdeApproximation.hpp"

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
SproutingRuleWithAnalyticalPdeApproximation<ELEMENT_DIM, SPACE_DIM>::SproutingRuleWithAnalyticalPdeApproximation(double MaxSproutingRate)
    : SproutingRule<ELEMENT_DIM, SPACE_DIM>(MaxSproutingRate)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SproutingRuleWithAnalyticalPdeApproximation<ELEMENT_DIM, SPACE_DIM>::~SproutingRuleWithAnalyticalPdeApproximation()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double SproutingRuleWithAnalyticalPdeApproximation<ELEMENT_DIM, SPACE_DIM>::GetSproutingProbability(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pParentCell){
    // initialisation 
    double Psprout;

    

    Psprout = mMaxSproutingRate*pow(vegf_concentration,1); // test since the concentration is between 0 and 1 

    return Psprout;
}

// Explicit instantiation
template class SproutingRuleWithAnalyticalPdeApproximation<1,1>;
template class SproutingRuleWithAnalyticalPdeApproximation<1,2>;
template class SproutingRuleWithAnalyticalPdeApproximation<1,3>;
template class SproutingRuleWithAnalyticalPdeApproximation<2,2>;
template class SproutingRuleWithAnalyticalPdeApproximation<2,3>;
template class SproutingRuleWithAnalyticalPdeApproximation<3,3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SproutingRuleWithAnalyticalPdeApproximation)