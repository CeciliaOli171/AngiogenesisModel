// #include "BloodFlowAngiogenesis.hpp"
// #include "ReplicatableVector.hpp"

// #include "Warnings.hpp"
// #include "Debug.hpp"

// template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// BloodFlowAngiogenesis<ELEMENT_DIM, SPACE_DIM>::BloodFlowAngiogenesis()
// {
// }

// template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// BloodFlowAngiogenesis<ELEMENT_DIM, SPACE_DIM>::~BloodFlowAngiogenesis()
// {
// }

// template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// void BloodFlowAngiogenesis<ELEMENT_DIM, SPACE_DIM>::SetFluxAtBoundaryNode(const Node<SPACE_DIM>& rNode, double flux)
// {
//     // check if the node is a boundary node 
//     // find the flux index 
//     // find the solution in the matrix / linear system
// }

// template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// void BloodFlowAngiogenesis<ELEMENT_DIM, SPACE_DIM>::SetPressureAtBoundaryNode(const Node<SPACE_DIM>& rNode, double pressure)
// {
//     // check if the node is a boundary node 
//     // find the pressure index 
//     // find the solution in the matrix / linear system
// }

// template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// double BloodFlowAngiogenesis<ELEMENT_DIM, SPACE_DIM>::GetFluxAtOutflow()
// {
//     // find outflow in vector solution
// }

// template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// void BloodFlowAngiogenesis<ELEMENT_DIM, SPACE_DIM>::GetSolutionAsFluxesAndPressures(std::vector<double>& rFluxesOnEdges, std::vector<double>& rPressuresOnNodes)
// {
// }

// template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// std::vector<double> BloodFlowAngiogenesis<ELEMENT_DIM, SPACE_DIM>::GetSolutionAsWSS(std::vector<double>& rFluxesOnEdges, std::vector<double>& rPressuresOnNodes)
// {
//     // find pressure, resistance and flow in solution vector 
// }

// template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// double BloodFlowAngiogenesis<ELEMENT_DIM, SPACE_DIM>::GetWSSAtNode(const Node<SPACE_DIM>& rNode)
// {
//     // calculates the WSS for one node 
// }

// // Explicit instantiation
// template class BloodFlowAngiogenesis<1,1>;
// template class BloodFlowAngiogenesis<1,2>;
// template class BloodFlowAngiogenesis<2,2>;
// template class BloodFlowAngiogenesis<1,3>;
// template class BloodFlowAngiogenesis<2,3>;
// template class BloodFlowAngiogenesis<3,3>;