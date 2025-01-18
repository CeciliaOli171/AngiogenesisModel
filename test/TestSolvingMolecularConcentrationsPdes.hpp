#ifndef TESTSOLVINGMOLECULARCONCENTRATIONSPDES_HPP_
#define TESTSOLVINGMOLECULARCONCENTRATIONSPDES_HPP_

// headliners
#include <cxxtest/TestSuite.h>
#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include "UblasIncludes.hpp"
#include "OutputFileHandler.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "Debug.hpp"
#include "CommandLineArguments.hpp"

// Geometry
#include "TetrahedralMesh.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "CellsGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "OffLatticeSimulation.hpp"
#include "AbstractCellPopulationBoundaryCondition.hpp"
#include "ChasteCuboid.hpp"

#include "PottsMesh.hpp"
#include "PottsMeshGenerator.hpp"
#include "PottsElement.hpp"
#include "PottsMeshReader.hpp"
#include "TrianglesMeshReader.hpp"

// Cell properties
#include "CellPropertyRegistry.hpp"
#include "CellLabel.hpp"
#include "CellData.hpp"
#include "CellId.hpp"
#include "Cell.hpp"
#include "CellBasedEventHandler.hpp"
#include "ApoptoticCellProperty.hpp"

#include "DifferentiatedCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "BranchingCellMutationState.hpp"
#include "TipCellMutationState.hpp"
#include "VesselCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "AbstractCellMutationState.hpp"

#include "FixedG1GenerationalCellCycleModel.hpp"
#include "UniformCellCycleModel.hpp"
#include "RandomDirectionCentreBasedDivisionRule.hpp"
#include "SproutingRule.hpp"
#include "DaughterCellModifier.hpp"

#include "PinnedCellsBoundaryCondition.hpp"

// Writers
#include "CellMutationStatesCountWriter.hpp"
#include "ConsecutiveBranchesWriter.hpp"
#include "BranchesNumberWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "BirthTimeCellWriter.hpp"
#include "TortuosityWriter.hpp"
#include "PottsMeshWriter.hpp"
#include "VtkMeshWriter.hpp"

// PDE solvers
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "SimpleNonlinearEllipticSolver.hpp"
#include "FunctionalBoundaryCondition.hpp"
#include "SimpleNewtonNonlinearSolver.hpp"
#include "VegfEquationPde.hpp"
#include "TrianglesMeshReader.hpp"
#include "EllipticBoxDomainPdeModifier.hpp"
#include "ParabolicBoxDomainPdeModifier.hpp"
#include "AveragedSourceEllipticPde.hpp"
#include "AveragedSourceParabolicPde.hpp"
#include "UniformSourceParabolicPde.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "MolecularConcentrationsDomainPdeModifier.hpp"
#include "CellwiseSourceParabolicPde.hpp"
#include "VegfBoundaryCondition.hpp"


class TestSolvingMolecularConcentrationsPdes : public AbstractCellBasedWithTimingsTestSuite
{
public:

    void NoTestSolvingVgefConcentrationPdeIn2D()
    {
        // Parameters input 
        CommandLineArguments* command_line = CommandLineArguments::Instance();
        double input_val_dudtcoeff = command_line->GetDoubleCorrespondingToOption("-dudtcoeff"); //1.0
        double input_val_diffusioncoeff = command_line->GetDoubleCorrespondingToOption("-diffusioncoeff"); //0.2
        double input_val_sourcecoeff = command_line->GetDoubleCorrespondingToOption("-sourcecoeff"); //1.0
        double input_val_creationcoeff = command_line->GetDoubleCorrespondingToOption("-creationcoeff"); //0.0
        double input_val_consumptioncoeff = command_line->GetDoubleCorrespondingToOption("-consumptioncoeff"); //0.0
        double input_val_initialvalue = command_line->GetDoubleCorrespondingToOption("-initialvalue"); //0.1
        double input_val_boundaryvalue = command_line->GetDoubleCorrespondingToOption("-boundaryvalue"); //0.1
        unsigned int input_val_endtime = command_line->GetDoubleCorrespondingToOption("-endtime"); 

        double boundary_cuboid_min = 0.0;
        double boundary_cuboid_max = 50.0;

        // creation of the mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0u, false, 27.0, 25.0));
        nodes.push_back(new Node<2>(1u, false, 26.0, 25.0));
        nodes.push_back(new Node<2>(2u, false, 25.0, 25.0));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingCellMutationState, p_branching_state); 
        MAKE_PTR(TipCellMutationState, p_tip_state);
        MAKE_PTR(VesselCellMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 
        MAKE_PTR(TransitCellProliferativeType, p_transit_type); // vessel segment ? 

        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        cell_population.Update(); // addition of this line compared to the sprouting test
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();

        // Set up simulation time for file output
        OffLatticeSimulation<2> simulator(cell_population);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, input_val_endtime);

        // Create PDE and boundary condition objects
        typedef VegfEquationPde<2> VegfEquationPde; 
        typedef VegfBoundaryCondition<2> VegfBoundaryCondition;
        MAKE_PTR_ARGS(VegfEquationPde, p_pde, (cell_population, input_val_dudtcoeff, input_val_diffusioncoeff, input_val_sourcecoeff, input_val_creationcoeff, input_val_consumptioncoeff));
        MAKE_PTR_ARGS(VegfBoundaryCondition, p_bc, (input_val_boundaryvalue, boundary_cuboid_min));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(boundary_cuboid_min, boundary_cuboid_min);
        ChastePoint<2> upper(boundary_cuboid_max, boundary_cuboid_max);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Initial conditions
        //std::vector<double> initial_condition_vector(2*mesh.GetNumNodes());
        // // tip cell has a vegf concentration of 0.1 for t=0
        // initial_condition_vector[4] = 1.0;
        // initial_condition_vector[5] = 1.0;
        //Vec initial_condition = PetscTools::CreateVec(initial_condition_vector);
        Vec initial_condition = nullptr;

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        typedef MolecularConcentrationsDomainPdeModifier<2> MolecularConcentrationsDomainPdeModifier;
        MAKE_PTR_ARGS(MolecularConcentrationsDomainPdeModifier, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 1.0, initial_condition, boundary_cuboid_min, input_val_initialvalue));
        p_pde_modifier->SetDependentVariableName("vegf_femesh_variable");
        p_pde_modifier->SetOutputGradient(true);

        p_pde_modifier->SetupSolve(cell_population,"TestVegfConcentrationPdeIn2D");

        // Run for 10 timesteps
        for (unsigned i=0; i<input_val_endtime; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
            p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
        }

        SimulationTime::Destroy();
    }


    void TestSolvingVgefConcentrationPdeIn3D()
    {
        // Parameters input 
        CommandLineArguments* command_line = CommandLineArguments::Instance();
        double input_val_dudtcoeff = command_line->GetDoubleCorrespondingToOption("-dudtcoeff"); //1.0
        double input_val_diffusioncoeff = command_line->GetDoubleCorrespondingToOption("-diffusioncoeff"); //0.2
        double input_val_sourcecoeff = command_line->GetDoubleCorrespondingToOption("-sourcecoeff"); //1.0
        double input_val_creationcoeff = command_line->GetDoubleCorrespondingToOption("-creationcoeff"); //0.0
        double input_val_consumptioncoeff = command_line->GetDoubleCorrespondingToOption("-consumptioncoeff"); //0.0
        double input_val_initialvalue = command_line->GetDoubleCorrespondingToOption("-initialvalue"); //0.1
        double input_val_boundaryvalue = command_line->GetDoubleCorrespondingToOption("-boundaryvalue"); //0.1
        unsigned int input_val_endtime = command_line->GetDoubleCorrespondingToOption("-endtime"); 

        double boundary_cuboid_min = 0.0;
        double boundary_cuboid_max = 50.0;

        // creation of the mesh
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u, false, 27.0, 25.0, 25.0));
        nodes.push_back(new Node<3>(1u, false, 26.0, 25.0, 25.0));
        nodes.push_back(new Node<3>(2u, false, 25.0, 25.0, 25.0));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingCellMutationState, p_branching_state); 
        MAKE_PTR(TipCellMutationState, p_tip_state);
        MAKE_PTR(VesselCellMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 
        MAKE_PTR(TransitCellProliferativeType, p_transit_type); // vessel segment ? 

        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        cell_population.Update(); // addition of this line compared to the sprouting test
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();

        // Set up simulation time for file output
        OffLatticeSimulation<3> simulator(cell_population);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(0.008333, input_val_endtime);

        // Create PDE and boundary condition objects
        typedef VegfEquationPde<3> VegfEquationPde; 
        typedef VegfBoundaryCondition<3> VegfBoundaryCondition;
        MAKE_PTR_ARGS(VegfEquationPde, p_pde, (cell_population, input_val_dudtcoeff, input_val_diffusioncoeff, input_val_sourcecoeff, input_val_creationcoeff, input_val_consumptioncoeff));
        MAKE_PTR_ARGS(VegfBoundaryCondition, p_bc, (input_val_boundaryvalue, boundary_cuboid_min));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<3> lower(boundary_cuboid_min, boundary_cuboid_min, boundary_cuboid_min);
        ChastePoint<3> upper(boundary_cuboid_max, boundary_cuboid_max, boundary_cuboid_max);
        MAKE_PTR_ARGS(ChasteCuboid<3>, p_cuboid, (lower, upper));

        // Initial conditions
        // std::vector<double> initial_condition_vector(3*mesh.GetNumNodes());
        // // tip cell has a vegf concentration of 0.1 for t=0
        // initial_condition_vector[6] = 0.1;
        // initial_condition_vector[7] = 0.1;
        // initial_condition_vector[8] = 0.1;
        // Vec initial_condition = PetscTools::CreateVec(initial_condition_vector);
        Vec initial_condition = nullptr;

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        typedef MolecularConcentrationsDomainPdeModifier<3> MolecularConcentrationsDomainPdeModifier;
        MAKE_PTR_ARGS(MolecularConcentrationsDomainPdeModifier, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 1.0, initial_condition, boundary_cuboid_min, input_val_initialvalue));
        p_pde_modifier->SetDependentVariableName("vegf_femesh_variable");
        p_pde_modifier->SetOutputGradient(true);

        p_pde_modifier->SetupSolve(cell_population,"TestVegfConcentrationPdeIn3D");

        // Run for 10 timesteps
        for (unsigned i=0; i<input_val_endtime; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
            p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
        }

        SimulationTime::Destroy();
    }
};

#endif /*TESTSOLVINGNONLINEARPDESTUTORIAL_HPP_*/
