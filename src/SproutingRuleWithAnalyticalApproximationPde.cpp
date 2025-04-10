#include "SproutingRuleWithAnalyticalApproximationPde.hpp"

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
SproutingRuleWithAnalyticalApproximationPde<ELEMENT_DIM, SPACE_DIM>::SproutingRuleWithAnalyticalApproximationPde(double MaxSproutingRateAnalyticalApproxPde, double diffusionCoefficient, double decayCoefficient, double creationCoefficient, double consumptionCoefficient, double boundaryCuboidMax, double sourceValue, double constantBackground, int PsproutFunctionTestNb)
    : SproutingRule<ELEMENT_DIM, SPACE_DIM>(MaxSproutingRateAnalyticalApproxPde), mMaxSproutingRateAnalyticalApproxPde(MaxSproutingRateAnalyticalApproxPde), mDiffusionCoefficient(diffusionCoefficient), mDecayCoefficient(decayCoefficient), mCreationCoefficient(creationCoefficient), mConsumptionCoefficient(consumptionCoefficient), mBoundaryCuboidMax(boundaryCuboidMax), mSourceValue(sourceValue), mConstantBackground(constantBackground), mPsproutFunctionTestNb(PsproutFunctionTestNb)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SproutingRuleWithAnalyticalApproximationPde<ELEMENT_DIM, SPACE_DIM>::~SproutingRuleWithAnalyticalApproximationPde()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double SproutingRuleWithAnalyticalApproximationPde<ELEMENT_DIM, SPACE_DIM>::GetVegfConcentrationAtNode(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pParentCell){
    c_vector<double, SPACE_DIM> x_parent = rCellPopulation.GetLocationOfCellCentre(pParentCell);

    double vegf_concentration = mSourceValue*(exp(-sqrt((mDecayCoefficient-mCreationCoefficient)/mDiffusionCoefficient)*x_parent[0]) - exp(sqrt((mDecayCoefficient-mCreationCoefficient)/mDiffusionCoefficient)*(x_parent[0]-2*mBoundaryCuboidMax)))/(1-exp(-sqrt(2*(mDecayCoefficient-mCreationCoefficient)/mDiffusionCoefficient)*mBoundaryCuboidMax))+mConstantBackground;

    return vegf_concentration;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double SproutingRuleWithAnalyticalApproximationPde<ELEMENT_DIM, SPACE_DIM>::GetSproutingProbability(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pParentCell){
    // initialisation 
    double Psprout;

    double vegf_concentration = GetVegfConcentrationAtNode(rCellPopulation, pParentCell);

    if(mPsproutFunctionTestNb == 0){
        // Linear function
        Psprout = mMaxSproutingRateAnalyticalApproxPde*vegf_concentration; // test since the concentration is between 0 and 1
    } else if(mPsproutFunctionTestNb == 1){
        // Hill function 
        double cmax = 1;
        double cmin = 0.5;
        double Pmax = 0.98;
        double Pmin = 0.4;
        double n = (1/log(cmax/cmin))*log((Pmax/Pmin)*(1-Pmin)/(1-Pmax));
        double K = pow(cmax*((1-Pmax)/Pmax),(1/n));
        Psprout = mMaxSproutingRateAnalyticalApproxPde*pow(vegf_concentration,n)/(pow(K, n) + pow(vegf_concentration,n));
    } 
     
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