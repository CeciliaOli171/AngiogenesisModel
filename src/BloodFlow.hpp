// #ifndef BLOODFLOW_HPP_
// #define BLOODFLOW_HPP_

// #include "AbstractCellPopulation.hpp"
// #include "TimeStepper.hpp"
// #include "VtkMeshWriter.hpp"

// /**
//  * A class for solving one-dimensional flow in pipe problems on branching trees.
//  *
//  * At graph edges: each pipe models Poiseuille flow (flux is linearly proportional to the pressure drop)
//  * At graph nodes: the flux is balanced so that mass is conserved.
//  *
//  * TO MODIFY!!!!!
//  * Works in 3D <1,3>
//  * Current functionality: pressure boundary conditions are set on each of the boundary nodes
//  * Solves for pressure at internal nodes and flux on edges
//  *
//  * In subclass: 
//  * * application to angiogenesis : all node pressures and edge fluxes are solved simultaneously using a direct matrix solution
//  */

// // use UpdateAtEndOfOutputTimeStep() or UpdateAtEndOfSolve() in DaughterCellModifier for blood flow?

// template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// class BloodFlow 
// {
// private:
//     // cell population -> go through it in Initialise to get a vascular network simplified for the calculation of blood flow etc 
//     AbstractCellPopulation<DIM> mrCellPopulation;
//     // outlet node i.e. parent vessel i.e. pinned cell location (have 2 branches coming from parent vessel and take one of them)
//     unsigned mOutletNodeIndex;
//     unsigned mInletNodeIndex;
//     // viscosity
//     double mBloodViscosity;
//     // boolean to choose if resistance changing dynamically or not 
//     bool mDynamicResistance; // false by default

//     double GetDynamicResistance();
//     double GetBloodViscosity();

//     // virtual void Initialise() (put it in public instead of private?): add in DaughterCellModifier that the branching cells have the list of the branches they are connected to; go through cell population and CellData to collect the branches, their length, the branching cells, the radius etc 
//     virtual void Initialise();

//     // add archive 

// public:
//     // constructor: inspired by constructor AbstractVentilationProblem, find boundarynode/give boundarynode + Initialise()
//     BloodFlow(AbstractCellPopulation<DIM>& rCellPopulation, unsigned outletNodeIndex, unsigned inletNodeIndex);

//     // destructor
//     virtual ~BloodFlow();

//     // get and set viscosity
//     double SetBloodViscosity();

//     // get and set resistance (virtual maybe to add in private) + (here virtual, in angiogenesis subclass constant resistance, add new subclass for varying resistance): based on radius, length and viscosity
//     double SetDynamicResistance();
//     double CalculateResistance(Element<ELEMENT_DIM,SPACE_DIM>& rElement, bool usePedley, double flux); // to modify variables and name

//     // calculates flux and pressure at boundary nodes (i.e. inlet, outlet and vessel tips)
//     virtual void SetFluxAtBoundaryNode(const Node<SPACE_DIM>& rNode, double flux)=0;
//     virtual void SetPressureAtBoundaryNode(const Node<SPACE_DIM>& rNode, double pressure)=0;

//     // calculates flux and pressure at vessel tip (important for vessel recession or growth)
//     void SetFluxAtVesselTipNode(const Node<SPACE_DIM>& rNode, double flux);
//     void SetPressureAtVesselTipNode(const Node<SPACE_DIM>& rNode, double pressure);

//     // Inflow
//     // SetConstantInflowPressures function to set the pressure at each inflow/leaf of the tree i.e. parent vessel in our case
//     // SetConstantInflowFluxes function to set the flux at each inflow/leaf of the tree i.e. parent vessel in our case
//     void SetConstantInflowFlux(double flux);
//     void SetConstantInflowPressure(double pressure);

//     // Outflow
//     // SetOutflowPressure function for pressure at parent vessel (virtual) related to SetConstantInFlowPressure
//     // SetOutflowFlux for flux at parent vessel (virtual) related to SetConstantInFlowFlux
//     void SetConstantOutflowFlux(double flux);
//     void SetConstantOutflowPressure(double pressure);

//     // function inspired by GetSolutionAsFluxesAndPressures to get fluxes and pressures ordered by edge index (virtual)
//     virtual void GetSolutionAsFluxesAndPressures(std::vector<double>& rFluxesOnEdges, std::vector<double>& rPressuresOnNodes)=0;

//     // calculates wall shear stresses and stores the results in a vector (virtual)
//     virtual std::vector<double> GetSolutionAsWSS(std::vector<double>& rFluxesOnEdges, std::vector<double>& rPressuresOnNodes)=0;

//     // access wall shear stress at a specific node? (from vector mSolution) + use pCell? (virtual)
//     virtual double GetWSSAtNode(const Node<SPACE_DIM>& rNode)=0;

//     // function solving the linear system inspired by SolveOverTime (check what flux balance at nodes and Poiseuille flow in the edges mean)
//     // may need to change names 
//     void SolveOverTime(TimeStepper& rTimeStepper, void (*pBoundaryConditionFunction)(AbstractVentilationProblem*, TimeStepper& rTimeStepper, const Node<3>&), const std::string& rDirName, const std::string& rFileBaseName);
    
//     // Solve function 
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