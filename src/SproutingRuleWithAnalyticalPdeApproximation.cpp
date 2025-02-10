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
SproutingRuleWithAnalyticalPdeApproximation<ELEMENT_DIM, SPACE_DIM>::SproutingRuleWithAnalyticalPdeApproximation(double MaxSproutingRateAnalyticalPde, double diffusionCoefficient, double decayCoefficient, double creationCoefficient, double consumptionCoefficient, double boundaryCuboidMax, double sourceValue, int testNb)
    : SproutingRule<ELEMENT_DIM, SPACE_DIM>(MaxSproutingRateAnalyticalPde), mMaxSproutingRateAnalyticalPde(MaxSproutingRateAnalyticalPde), mDiffusionCoefficient(diffusionCoefficient), mDecayCoefficient(decayCoefficient), mCreationCoefficient(creationCoefficient), mConsumptionCoefficient(consumptionCoefficient), mBoundaryCuboidMax(boundaryCuboidMax), mSourceValue(sourceValue), mTestNb(testNb)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SproutingRuleWithAnalyticalPdeApproximation<ELEMENT_DIM, SPACE_DIM>::~SproutingRuleWithAnalyticalPdeApproximation()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double SproutingRuleWithAnalyticalPdeApproximation<ELEMENT_DIM, SPACE_DIM>::GetVegfConcentrationAtNode(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pParentCell){
    c_vector<double, SPACE_DIM> x_parent = rCellPopulation.GetLocationOfCellCentre(pParentCell);

    double vegf_concentration = mSourceValue*(exp(-sqrt((mDecayCoefficient-mCreationCoefficient)/mDiffusionCoefficient)*x_parent[0]) - exp(sqrt((mDecayCoefficient-mCreationCoefficient)/mDiffusionCoefficient)*(x_parent[0]-2*mBoundaryCuboidMax)))/(1-exp(-sqrt(2*(mDecayCoefficient-mCreationCoefficient)/mDiffusionCoefficient)*mBoundaryCuboidMax));

    return vegf_concentration;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double SproutingRuleWithAnalyticalPdeApproximation<ELEMENT_DIM, SPACE_DIM>::GetSproutingProbability(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pParentCell){
    // initialisation 
    double Psprout;

    double vegf_concentration = GetVegfConcentrationAtNode(rCellPopulation, pParentCell);

    if(mTestNb == 0){
        // Linear function
        Psprout = mMaxSproutingRateAnalyticalPde*vegf_concentration; // test since the concentration is between 0 and 1
    } else if(mTestNb == 1){
        // Hill function (n=1) : Witzenbichler et al. 1998
        Psprout = mMaxSproutingRateAnalyticalPde*pow(vegf_concentration,1)/(pow(0.23, 1) + pow(vegf_concentration,1));
    } else if(mTestNb == 2){
        // Hill function (n=2) : Witzenbichler et al. 1998
        Psprout = mMaxSproutingRateAnalyticalPde*pow(vegf_concentration,2)/(pow(0.23, 2) + pow(vegf_concentration,2));
    } else if(mTestNb == 3){
        // Hill function (n=3) : Witzenbichler et al. 1998
        Psprout = mMaxSproutingRateAnalyticalPde*pow(vegf_concentration,3)/(pow(0.23, 3) + pow(vegf_concentration,3));
    }
     
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