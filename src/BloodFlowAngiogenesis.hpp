// #ifndef BLOODFLOWANGIOGENESIS_HPP_
// #define BLOODFLOWANGIOGENESIS_HPP_

// #include "TetrahedralMesh.hpp"
// #include "TimeStepper.hpp"
// #include "VtkMeshWriter.hpp"
// #include "BloodFlow.hpp"

// /**
//  * TO MODIFY!!!! 
//  * A class for solving one-dimensional flow in pipe problems on branching trees.
//  *
//  * At graph edges: each pipe models Poiseuille flow (flux is linearly proportional to the pressure drop)
//  * At graph nodes: the flux is balanced so that mass is conserved.
//  *
//  * Works in 3D <1,3>
//  * Current functionality: pressure boundary conditions are set on each of the boundary nodes
//  * Solves for pressure at internal nodes and flux on edges
//  *
//  * In this subclass all node pressures and edge fluxes are solved simultaneously using a direct matrix solution.
//  */

// template class<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// class BloodFlowAngiogenesis : public BloodFlow 
// {
//     friend class TestSolvingBloodFlow;

// private:
//     // linear system 
//     LinearSystem* mpLinearSystem;
//     // mFluxScaling for pressure and flux solution at a comparable magnitude?
//     // vector solution 
//     Vec mSolution;

//     // Assemble() function with the option of considering a dynamic resistance (necessary for us since angiogenesis = branches lengths vary)
//     void Assemble(bool dynamicReassemble=false);

//     // add archive 

// public:
//     // constructor 
//     BloodFlowAngiogenesis();

//     // destructor: destroy linear system (important to not forget)
//     ~BloodFlowAngiogenesis();

//     // calculates flux and pressure at boundary nodes: use pCell instead of Node?
//     void SetFluxAtBoundaryNode(const Node<SPACE_DIM>& rNode, double flux);
//     void SetPressureAtBoundaryNode(const Node<SPACE_DIM>& rNode, double pressure);

//     // returns the flux for the outflow (already in the BloodFlow class)
//     double GetFluxAtOutflow();

//     // Solver
//     void Solve();

//     // organises the vector solution into flux and pressures components
//     void GetSolutionAsFluxesAndPressures(std::vector<double>& rFluxesOnEdges, std::vector<double>& rPressuresOnNodes);

//     // calculates wall shear stresses and stores the results in a vector
//     std::vector<double> GetSolutionAsWSS(std::vector<double>& rFluxesOnEdges, std::vector<double>& rPressuresOnNodes);

//     // access wall shear stress at a specific node? (from vector mSolution) + use pCell?
//     double GetWSSAtNode(const Node<SPACE_DIM>& rNode);
// };

// #endif /* BLOODFLOWANGIOGENESIS_HPP_ */