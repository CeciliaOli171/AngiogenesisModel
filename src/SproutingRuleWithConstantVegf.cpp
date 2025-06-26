#include "SproutingRuleWithConstantVegf.hpp"

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
#include <cmath>
#include "ReplicatableVector.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SproutingRuleWithConstantVegf<ELEMENT_DIM, SPACE_DIM>::SproutingRuleWithConstantVegf(double MaxSproutingRateConstantVegf, double diffusionCoefficient, double decayCoefficient, double creationCoefficient, double consumptionCoefficient, double sourceValue, double constantBackground, double cMax, double cMin, double pMax, double pMin, double boundaryCuboidMax, int PsproutFunctionTestNb)
    : SproutingRule<ELEMENT_DIM, SPACE_DIM>(MaxSproutingRateConstantVegf), mMaxSproutingRateConstantVegf(MaxSproutingRateConstantVegf), mConstantBackground(constantBackground)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SproutingRuleWithConstantVegf<ELEMENT_DIM, SPACE_DIM>::~SproutingRuleWithConstantVegf()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double SproutingRuleWithConstantVegf<ELEMENT_DIM, SPACE_DIM>::GetVegfConcentrationAtNode(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pParentCell){
    return mConstantBackground;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double SproutingRuleWithConstantVegf<ELEMENT_DIM, SPACE_DIM>::GetSproutingProbability(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pParentCell){
    // initialisation 
    double Psprout;

    double vegf_concentration = GetVegfConcentrationAtNode(rCellPopulation, pParentCell);

    if(mPsproutFunctionTestNb == 0){
        // Linear function
        Psprout = mMaxSproutingRateAnalyticalApproxPde*vegf_concentration; // test since the concentration is between 0 and 1
    } else if(mPsproutFunctionTestNb == 1){
        // Hill function 
        double n = (1/log(mCMax/mCMin))*log((mPMax/mPMin)*(1-mPMin)/(1-mPMax));
        double K = pow(mCMax*((1-mPMax)/mPMax),(1/n));
        Psprout = mMaxSproutingRateAnalyticalApproxPde*pow(vegf_concentration,n)/(pow(K, n) + pow(vegf_concentration,n));
    } 

    mMaxSproutingRateConstantVegf = Psprout;
     
    return mMaxSproutingRateConstantVegf;
}

// Explicit instantiation
template class SproutingRuleWithConstantVegf<1,1>;
template class SproutingRuleWithConstantVegf<1,2>;
template class SproutingRuleWithConstantVegf<1,3>;
template class SproutingRuleWithConstantVegf<2,2>;
template class SproutingRuleWithConstantVegf<2,3>;
template class SproutingRuleWithConstantVegf<3,3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SproutingRuleWithConstantVegf)