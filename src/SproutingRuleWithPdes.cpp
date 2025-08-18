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
#include <cmath>
#include "ReplicatableVector.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SproutingRuleWithPdes<ELEMENT_DIM, SPACE_DIM>::SproutingRuleWithPdes(double MaxSproutingRatePdes, double thresholdLength, boost::shared_ptr<AbstractBoxDomainPdeModifier<SPACE_DIM> > pPdeModifier, int PsproutFunctionTestNb)
    : SproutingRule<ELEMENT_DIM, SPACE_DIM>(MaxSproutingRatePdes, thresholdLength), mMaxSproutingRatePdes(MaxSproutingRatePdes), mpPdeModifier(pPdeModifier), mPsproutFunctionTestNb(PsproutFunctionTestNb)
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

    if(mPsproutFunctionTestNb == 0){
        // Linear function
        Psprout = mMaxSproutingRatePdes*vegf_concentration; // test since the concentration is between 0 and 1
    } else if(mPsproutFunctionTestNb == 1){
        // Hill function 
        double cmax = 1;
        double cmin = 0.5;
        double Pmax = 0.98;
        double Pmin = 0.4;
        double n = (1/log(cmax/cmin))*log((Pmax/Pmin)*(1-Pmin)/(1-Pmax));
        double K = pow(cmax*((1-Pmax)/Pmax),(1/n));
        Psprout = mMaxSproutingRatePdes*pow(vegf_concentration,n)/(pow(K, n) + pow(vegf_concentration,n));
    }

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