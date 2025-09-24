#include "SproutingRuleTwoLesionsYAxis.hpp"

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
SproutingRuleTwoLesionsYAxis<ELEMENT_DIM, SPACE_DIM>::SproutingRuleTwoLesionsYAxis(double MaxSproutingRateTwoLesionsYAxis, double thresholdLength, double diffusionCoefficient, double decayCoefficient, double creationCoefficient, double consumptionCoefficient, double sourceValue1, double sourceValue2, double constantBackground, double cMax, double cMin, double pMax, double pMin, double locationLesion1, double locationLesion2, int PsproutFunctionTestNb)
    : SproutingRule<ELEMENT_DIM, SPACE_DIM>(MaxSproutingRateTwoLesionsYAxis, thresholdLength), mMaxSproutingRateTwoLesionsYAxis(MaxSproutingRateTwoLesionsYAxis), mDiffusionCoefficient(diffusionCoefficient), mDecayCoefficient(decayCoefficient), mCreationCoefficient(creationCoefficient), mConsumptionCoefficient(consumptionCoefficient), mSourceValue1(sourceValue1), mSourceValue2(sourceValue2), mConstantBackground(constantBackground), mCMax(cMax), mCMin(cMin), mPMax(pMax), mPMin(pMin), mLocationLesion1(locationLesion1), mLocationLesion2(locationLesion2), mPsproutFunctionTestNb(PsproutFunctionTestNb)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SproutingRuleTwoLesionsYAxis<ELEMENT_DIM, SPACE_DIM>::~SproutingRuleTwoLesionsYAxis()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double SproutingRuleTwoLesionsYAxis<ELEMENT_DIM, SPACE_DIM>::GetVegfConcentrationAtNode(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pParentCell){
    c_vector<double, SPACE_DIM> x_parent = rCellPopulation.GetLocationOfCellCentre(pParentCell);

    double Kc = sqrt((mDecayCoefficient-mCreationCoefficient)/mDiffusionCoefficient);

    double vegf_concentration_1 = (mSourceValue1-mConstantBackground)*exp(-Kc*(x_parent[1]-mLocationLesion1))+mConstantBackground;
    double vegf_concentration_2 = (mSourceValue2-mConstantBackground)*exp(Kc*(x_parent[1]-mLocationLesion2))+mConstantBackground;

    double vegf_concentration = (vegf_concentration_1+vegf_concentration_2)/2;

    return vegf_concentration;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double SproutingRuleTwoLesionsYAxis<ELEMENT_DIM, SPACE_DIM>::GetSproutingProbability(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pParentCell){
    // initialisation 
    double Psprout;

    double vegf_concentration = GetVegfConcentrationAtNode(rCellPopulation, pParentCell);

    if(mPsproutFunctionTestNb == 0){
        // Linear function
        Psprout = mMaxSproutingRateTwoLesionsYAxis*vegf_concentration; // test since the concentration is between 0 and 1
    } else if(mPsproutFunctionTestNb == 1){
        // Hill function 
        double n = (1/log(mCMax/mCMin))*log((mPMax/mPMin)*(1-mPMin)/(1-mPMax));
        double K = pow(mCMax*((1-mPMax)/mPMax),(1/n));
        Psprout = mMaxSproutingRateTwoLesionsYAxis*pow(vegf_concentration,n)/(pow(K, n) + pow(vegf_concentration,n));
    } 
    
    return Psprout;
}

// Explicit instantiation
template class SproutingRuleTwoLesionsYAxis<1,1>;
template class SproutingRuleTwoLesionsYAxis<1,2>;
template class SproutingRuleTwoLesionsYAxis<1,3>;
template class SproutingRuleTwoLesionsYAxis<2,2>;
template class SproutingRuleTwoLesionsYAxis<2,3>;
template class SproutingRuleTwoLesionsYAxis<3,3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SproutingRuleTwoLesionsYAxis)