// #include "BloodFlowAngiogenesis.hpp"
// #include "ReplicatableVector.hpp"

// #include "Warnings.hpp"
// #include "Debug.hpp"

// template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// BloodFlowAngiogenesis<ELEMENT_DIM, SPACE_DIM>::BloodFlowAngiogenesis(): AbstractVentilationProblem(rMeshDirFilePath, rootIndex), mpLinearSystem(nullptr), mSolution(nullptr)
// {
//     // We solve for flux at every edge and for pressure at each node/bifurcation
//     // Note pipe flow equation has 3 variables and flux balance has 3 variables (at a bifurcation)
//     // preallocating 5 non-zeros allows for 4-way branching
//     mSolution = PetscTools::CreateVec(mMesh.GetNumNodes()+mMesh.GetNumElements());
//     mpLinearSystem = new LinearSystem(mSolution, 5u);
//     mpLinearSystem->SetAbsoluteTolerance(1e-5);
// }

// template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// BloodFlowAngiogenesis<ELEMENT_DIM, SPACE_DIM>::~BloodFlowAngiogenesis()
// {
//     if (mpLinearSystem)
//     {
//         delete mpLinearSystem;
//     }
//     if (mSolution)
//     {
//         PetscTools::Destroy(mSolution);
//     }
// }

// template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// void BloodFlowAngiogenesis<ELEMENT_DIM, SPACE_DIM>::Assemble(bool dynamicReassemble)
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

// template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// void BloodFlowAngiogenesis<ELEMENT_DIM, SPACE_DIM>::Solve()
// {
// }

// // Explicit instantiation
// template class BloodFlowAngiogenesis<1,1>;
// template class BloodFlowAngiogenesis<1,2>;
// template class BloodFlowAngiogenesis<2,2>;
// template class BloodFlowAngiogenesis<1,3>;
// template class BloodFlowAngiogenesis<2,3>;
// template class BloodFlowAngiogenesis<3,3>;