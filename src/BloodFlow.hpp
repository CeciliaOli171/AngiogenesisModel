// #ifndef BLOODFLOW_HPP_
// #define BLOODFLOW_HPP_

// #include "AbstractCellPopulation.hpp"
// #include "TimeStepper.hpp"
// #include "VtkMeshWriter.hpp"

// /**
//  * A class for solving one-dimensional flow in pipe problems on branching trees.
//  *
//  * graph = vascular tree 
//  * At graph edges: each pipe models Poiseuille flow (flux is linearly proportional to the pressure drop)
//  * At graph nodes = branching cells: the flux is balanced so that mass is conserved.
//  *
//  * TO MODIFY!!!!!
//  * Solves for pressure at internal nodes and flux on edges
//  *
//  * In subclass: 
//  * * application to angiogenesis : all node pressures and edge fluxes are solved simultaneously using a direct matrix solution
//  */

// // use UpdateAtEndOfOutputTimeStep() or UpdateAtEndOfSolve() in DaughterCellModifier for blood flow?
// // branching cells should have inside CellData the branch numbers of the branches it is connected to (to construct the matrix)

// template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// class BloodFlow 
// {
// private:
//     /* Parameters */
//     AbstractCellPopulation<DIM> mrCellPopulation; // cell population -> go through it in Initialise to get a vascular network simplified for the calculation of blood flow etc 
//     unsigned mOutletNodeIndex; // outlet node i.e. in parent vessel i.e. one of the pinned cell location 
//     unsigned mInletNodeIndex; // inlet node i.e. in parent vessel i.e. the other pinned cell location
//     double mBloodViscosity; // viscosity
//     bool mDynamicResistance; // boolean to choose if resistance changing dynamically or not, false by default 

//     /**
//      * @return the blood viscosity of the vascular network in kg/(mm*sec).
//      */
//     double GetBloodViscosity();

//     /**
//      * @return the boolean deciding if the resistance is dynamic (flux related) or static (constant).
//      */
//     bool GetDynamicResistance();

//     // virtual void Initialise() (put it in public instead of private?): add in DaughterCellModifier that the branching cells have the list of the branches they are connected to; go through cell population and CellData to collect the branches, their length, the branching cells, the radius etc 
//     virtual void Initialise();

//     // add archive ?

// public:
//     // constructor: inspired by constructor AbstractVentilationProblem, find boundarynode/give boundarynode + Initialise() by creating a linear system of the appropriate size to match the mesh
//     BloodFlow(AbstractCellPopulation<DIM>& rCellPopulation, unsigned outletNodeIndex, unsigned inletNodeIndex);

//     /**
//      * Destructor.
//      */
//     virtual ~BloodFlow();

//     /**
//      * Set the blood viscosity in the vascular network.
//      * 
//      * @param bloodViscosity the viscosity in kg/(mm*sec)
//      *
//      */
//     double SetBloodViscosity(double bloodViscosity);

//     /**
//      * Set if the resistance inside a branch is dynamic (flux dependent) or static (constant).
//      *
//      * @param dynamicResistance the boolean associated to the dynamic resistance.
//      */
//     double SetDynamicResistance(bool dynamicResistance);

//     // maybe virtual
//     // constant in this class and add an option for dynamic in angiogenesis class? 
//     double CalculateResistance(Element<ELEMENT_DIM,SPACE_DIM>& rElement, bool usePedley, double flux); // to modify variables and name

//     /**
//      * Sets a flux boundary condition for a given node (inlet, outlet or vessel tip).
//      *
//      * The given boundary condition will be applied at the next time step and persist through
//      * time unless overwritten.
//      *
//      * @param rNode The node to set the boundary condition for
//      * @param flux The flux boundary condition in (m^3)/s
//      */
//     virtual void SetFluxAtBoundaryNode(const Node<SPACE_DIM>& rNode, double flux)=0;

//     /**
//      * Sets a pressure boundary condition for a given node (inlet, outler or vessel tip).
//      *
//      * The given boundary condition will be applied at the next time step and persist through
//      * time unless overwritten.
//      *
//      * @param rNode The node to set the boundary condition for
//      * @param pressure The pressure boundary condition in Pascals
//      */
//     virtual void SetPressureAtBoundaryNode(const Node<SPACE_DIM>& rNode, double pressure)=0;

//     // calculates flux and pressure at vessel tip (important for vessel recession or growth)
//     void SetFluxAtVesselTipNode(const Node<SPACE_DIM>& rNode, double flux);
//     void SetPressureAtVesselTipNode(const Node<SPACE_DIM>& rNode, double pressure);

//     // Inflow
//     // SetConstantInflowPressures function to set the pressure at each inflow/leaf of the tree i.e. inlet parent vessel in our case
//     // SetConstantInflowFluxes function to set the flux at each inflow/leaf of the tree i.e. parent vessel in our case
//     void SetConstantInflowFlux(double flux);
//     void SetConstantInflowPressure(double pressure);

//     // Outflow
//     // SetOutflowPressure function for pressure at parent vessel (virtual) related to SetConstantInFlowPressure
//     // SetOutflowFlux for flux at parent vessel (virtual) related to SetConstantInFlowFlux
//     void SetConstantOutflowFlux(double flux);
//     void SetConstantOutflowPressure(double pressure);

//     /**
//      * Calculates the solution as a vector of fluxes and pressures ordered by edges and nodes respectively.
//      *
//      * @param rFluxesOnEdges The fluxes ordered by edge index (this vector is resized)
//      * @param rPressuresOnNodes The pressures ordered by node index (this vector is resized)
//      */
//     virtual void GetSolutionAsFluxesAndPressures(std::vector<double>& rFluxesOnEdges, std::vector<double>& rPressuresOnNodes)=0;

//     /**
//      * Calculates the wall shear stresses from the fluxes and pressures obtained in the solution vector.
//      *
//      * @param rFluxesOnEdges The fluxes ordered by edge index (this vector is resized)
//      * @param rPressuresOnNodes The pressures ordered by node index (this vector is resized)
//      */
//     virtual std::vector<double> GetSolutionAsWSS(std::vector<double>& rFluxesOnEdges, std::vector<double>& rPressuresOnNodes)=0;

//     // access wall shear stress at a specific node or a specific cell? (from vector mSolution) + use pCell? (virtual)
//     virtual double GetWSSAtNode(const Node<SPACE_DIM>& rNode)=0;

//     void Initialise();

//     // function solving the linear system inspired by SolveOverTime (check what flux balance at nodes and Poiseuille flow in the edges mean)
//     // may need to change names 
//     void SolveOverTime(TimeStepper& rTimeStepper, void (*pBoundaryConditionFunction)(AbstractVentilationProblem*, TimeStepper& rTimeStepper, const Node<3>&), const std::string& rDirName, const std::string& rFileBaseName);
    
//     /**
//      *  Solve the system at one time-point. 
//      */ 
//     virtual void Solve()=0;

// // writes a VTK file 
// #ifdef CHASTE_VTK
//     /**
//      * Add flux and pressure data to a VtkMeshWriter.
//      * @param rVtkWriter  the mesh writer ready for the data
//      * @param rSuffix  Suffix with which to annotate e.g. pressure_001
//      */
//     void AddDataToVtk(VtkMeshWriter<1, 3>& rVtkWriter, const std::string& rSuffix);

//     /**
//      * Output the solution to a Vtk file
//      * @param rDirName A directory name relative to CHASTE_TEST_OUTPUT.
//      * @param rFileBaseName The base name of the new VTK file.
//      */
//     void WriteVtk(const std::string& rDirName, const std::string& rFileBaseName);
// #endif // CHASTE_VTK

// };

// #endif /* BLOODFLOW_HPP_ */