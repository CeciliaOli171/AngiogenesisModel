#include "SproutingRuleWithAnalyticalApproximationPde.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "VesselSegmentMutationState.hpp"
#include "VesselTipMutationState.hpp"
#include "BranchingSegmentMutationState.hpp"
#include "UniformCellCycleModel.hpp"

#include "RandomNumberGenerator.hpp"
#include <algorithm>
#include "Debug.hpp"
#include <random>
#include <cstdlib>
#include <iterator>
#include <cmath>
#include "ReplicatableVector.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SproutingRuleWithAnalyticalApproximationPde<ELEMENT_DIM, SPACE_DIM>::SproutingRuleWithAnalyticalApproximationPde(double MaxSproutingRateAnalyticalApproxPde, double thresholdLength, double diffusionCoefficient, double decayCoefficient, double creationCoefficient, double consumptionCoefficient, double sourceValue, double constantBackground, double cMax, double cMin, double pMax, double pMin)
    : SproutingRule<ELEMENT_DIM, SPACE_DIM>(MaxSproutingRateAnalyticalApproxPde, thresholdLength), mMaxSproutingRateAnalyticalApproxPde(MaxSproutingRateAnalyticalApproxPde), mDiffusionCoefficient(diffusionCoefficient), mDecayCoefficient(decayCoefficient), mCreationCoefficient(creationCoefficient), mConsumptionCoefficient(consumptionCoefficient), mSourceValue(sourceValue), mConstantBackground(constantBackground), mCMax(cMax), mCMin(cMin), mPMax(pMax), mPMin(pMin)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SproutingRuleWithAnalyticalApproximationPde<ELEMENT_DIM, SPACE_DIM>::~SproutingRuleWithAnalyticalApproximationPde()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double SproutingRuleWithAnalyticalApproximationPde<ELEMENT_DIM, SPACE_DIM>::GetVegfConcentrationAtNode(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pParentCell){
    c_vector<double, SPACE_DIM> x_parent = rCellPopulation.GetLocationOfCellCentre(pParentCell);

    double Kc = sqrt((mDecayCoefficient-mCreationCoefficient)/mDiffusionCoefficient);
    double vegf_concentration = (mSourceValue-mConstantBackground)*exp(-Kc*abs(x_parent[0]))+mConstantBackground;

    return vegf_concentration;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double SproutingRuleWithAnalyticalApproximationPde<ELEMENT_DIM, SPACE_DIM>::GetSproutingProbability(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pParentCell){
    // initialisation 
    double Psprout;

    double vegf_concentration = GetVegfConcentrationAtNode(rCellPopulation, pParentCell);

    // Hill function 
    double n = (1/log(mCMax/mCMin))*log((mPMax/mPMin)*(1-mPMin)/(1-mPMax));
    double K = pow(mCMax*((1-mPMax)/mPMax),(1/n));
    Psprout = mMaxSproutingRateAnalyticalApproxPde*pow(vegf_concentration,n)/(pow(K, n) + pow(vegf_concentration,n));
     
    return Psprout;
}

// Explicit instantiation
template class SproutingRuleWithAnalyticalApproximationPde<1,1>;
template class SproutingRuleWithAnalyticalApproximationPde<1,2>;
template class SproutingRuleWithAnalyticalApproximationPde<1,3>;
template class SproutingRuleWithAnalyticalApproximationPde<2,2>;
template class SproutingRuleWithAnalyticalApproximationPde<2,3>;
template class SproutingRuleWithAnalyticalApproximationPde<3,3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SproutingRuleWithAnalyticalApproximationPde)