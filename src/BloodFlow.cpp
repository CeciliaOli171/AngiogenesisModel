// #include "BloodFlow.hpp"

// #include "VesselTipMutationState.hpp"

// template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// BloodFlow<ELEMENT_DIM, SPACE_DIM>::BloodFlow(AbstractCellPopulation<DIM>& rCellPopulation, unsigned outletNodeIndex, unsigned inletNodeIndex): mrCellPopulation(rCellPopulation), mOutletNodeIndex(outletNodeIndex), mInletNodeIndex(inletNodeIndex), mBloodViscosity(), mDynamicResistance(false)
// {
//     // add something to check that the outlet node is a boundary node?
//     Initialise();
// }

// template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// void BloodFlow<ELEMENT_DIM, SPACE_DIM>::Initialise()
// {
// }

// template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// double BloodFlow<ELEMENT_DIM, SPACE_DIM>::CalculateResistance(Element<1,3>& rElement, bool usePedley, double flux)
// {
// }

// template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// void BloodFlow<ELEMENT_DIM, SPACE_DIM>::SetFluxAtVesselTipNode(double flux)
// {
//     for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = mrCellPopulation.Begin(); cell_iter != mrCellPopulation.End(); ++cell_iter)
//     {
//         if ((*cell_iter)->GetMutationState()->template IsType<VesselTipMutationState>()) 
//         {
//             SetFluxAtBoundaryNode(*(*cell_iter), flux);
//         }
//     }
// }

// template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// void BloodFlow<ELEMENT_DIM, SPACE_DIM>::SetPressureAtVesselTipNode(double pressure)
// {
//     for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = mrCellPopulation.Begin(); cell_iter != mrCellPopulation.End(); ++cell_iter)
//     {
//         if ((*cell_iter)->GetMutationState()->template IsType<VesselTipMutationState>()) 
//         {
//             SetPressureAtBoundaryNode(*(*cell_iter), pressure);
//         }
//     }
// }

// template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// void BloodFlow<ELEMENT_DIM, SPACE_DIM>::SetConstantInflowFlux(double flux)
// {
//     SetFluxAtBoundaryNode(*(mrCellPopulation.GetNode(mInletNodeIndex)), flux);
// }

// template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// void BloodFlow<ELEMENT_DIM, SPACE_DIM>::SetConstantInflowPressure(double pressure)
// {
//     SetPressureAtBoundaryNode(*(mrCellPopulation.GetNode(mInletNodeIndex)), pressure);
// }

// template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// void BloodFlow<ELEMENT_DIM, SPACE_DIM>::SetConstantOutflowFlux(double flux)
// {
//     SetFluxAtBoundaryNode(*(mrCellPopulation.GetNode(mOutletNodeIndex)), flux);
// }

// template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// void BloodFlow<ELEMENT_DIM, SPACE_DIM>::SetConstantOutflowPressure(double pressure)
// {
//     SetPressureAtBoundaryNode(*(mrCellPopulation.GetNode(mOutletNodeIndex)), pressure);
// }

// template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// void BloodFlow<ELEMENT_DIM, SPACE_DIM>::SolveOverTime()
// {
//     // to complete
// }

// // Explicit instantiation
// template class BloodFlow<1,1>;
// template class BloodFlow<1,2>;
// template class BloodFlow<2,2>;
// template class BloodFlow<1,3>;
// template class BloodFlow<2,3>;
// template class BloodFlow<3,3>;