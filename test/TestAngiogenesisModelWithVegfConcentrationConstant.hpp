#ifndef TESTANGIOGENESISMODELWITHVEGFCONCENTRATIONCONSTANT_HPP_
#define TESTANGIOGENESISMODELWITHVEGFCONCENTRATIONCONSTANT_HPP_

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
#include "BranchingSegmentMutationState.hpp"
#include "VesselTipMutationState.hpp"
#include "VesselSegmentMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "AbstractCellMutationState.hpp"

#include "FixedG1GenerationalCellCycleModel.hpp"
#include "UniformCellCycleModel.hpp"
#include "RandomDirectionCentreBasedDivisionRule.hpp"
#include "SproutingRule.hpp"
#include "SproutingRuleWithConstantVegf.hpp"
#include "SproutingRuleWithPdes.hpp"
#include "SproutingRuleWithAnalyticalApproximationPde.hpp"
#include "DaughterCellModifier.hpp"

#include "PinnedCellsBoundaryCondition.hpp"

// Forces 
#include "GeneralisedLinearSpringForce.hpp"
#include "DifferentialAdhesionGeneralisedLinearSpringForce.hpp"
#include "ChemoForce.hpp"
#include "RandomForce.hpp"
#include "PersistenceForce.hpp"
#include "AngularForce.hpp"
#include "LinearMechanicalForceModified.hpp"
#include "DirectionalPersistenceCellModifier.hpp"
#include "ChemoForceWithPdes.hpp"
#include "ChemoForceWithAnalyticalApproximationPde.hpp"
#include "ChemoForceWithConstantVegf.hpp"

// Writers
#include "CellMutationStatesCountWriter.hpp"
#include "ConsecutiveBranchesWriter.hpp"
#include "BranchNumberWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "BirthTimeCellWriter.hpp"
#include "TortuosityWriter.hpp"
#include "PottsMeshWriter.hpp"
#include "VtkMeshWriter.hpp"
#include "AnastomosisWriter.hpp"
#include "ConnectivityWriter.hpp"

// PDE solvers
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "SimpleNonlinearEllipticSolver.hpp"
#include "FunctionalBoundaryCondition.hpp"
#include "SimpleNewtonNonlinearSolver.hpp"
#include "VegfEquationPde.hpp"
#include "EllipticBoxDomainPdeModifier.hpp"
#include "ParabolicBoxDomainPdeModifier.hpp"
#include "AveragedSourceEllipticPde.hpp"
#include "AveragedSourceParabolicPde.hpp"
#include "UniformSourceParabolicPde.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "MolecularConcentrationsBoxDomainPdeModifier.hpp"
#include "CellwiseSourceParabolicPde.hpp"
#include "VegfBoundaryCondition.hpp"



class TestAngiogenesisModelWithVegfConcentrationConstant : public AbstractCellBasedTestSuite
{

public:

    void TestAngiogenesisModelWithVegfConcentrationConstantIn2D() 
    {
        ///////////
        // INPUT // 
        ///////////

        // Parameters input 
        CommandLineArguments* command_line = CommandLineArguments::Instance();

        // parameters for vegf concentration
        double input_val_vegf_constantbackground = command_line->GetDoubleCorrespondingToOption("-vegfconstantbackground"); //0.1

        // parameters for forces 
        double input_val_S = command_line->GetDoubleCorrespondingToOption("-S"); //1E-5
        double input_val_sigma = command_line->GetDoubleCorrespondingToOption("-sigma"); //1E-1
        double input_val_chi = command_line->GetDoubleCorrespondingToOption("-chi"); //1E-4
        double input_val_omegap = command_line->GetDoubleCorrespondingToOption("-omegap"); //1E-1
        double input_val_omegaa = command_line->GetDoubleCorrespondingToOption("-omegaa"); //1E-5

        // parameters for Psprout 
        double input_val_maxsproutingrate = command_line->GetDoubleCorrespondingToOption("-maxsproutingrate"); 
        double input_val_cmax = command_line->GetDoubleCorrespondingToOption("-cmax");
        double input_val_cmin = command_line->GetDoubleCorrespondingToOption("-cmin");
        double input_val_pmax = command_line->GetDoubleCorrespondingToOption("-pmax");
        double input_val_pmin = command_line->GetDoubleCorrespondingToOption("-pmin");

        // parameters for anastomosis 
        double input_val_anastomosislength = command_line->GetDoubleCorrespondingToOption("-anastomosislength");
        double input_val_thresholdlength = command_line->GetDoubleCorrespondingToOption("-thresholdlength");;

        // general parameters (time, random seed, output directory)
        double input_val_time = command_line->GetDoubleCorrespondingToOption("-time"); 
        double input_val_seed = command_line->GetIntCorrespondingToOption("-seed"); 
        std::string output_directory_vegfconstant = command_line->GetStringCorrespondingToOption("-output_directory_vegfconstant");

        ////////////
        // MESHES // 
        ////////////

        // mesh for pdes 
        //double boundary_cuboid_min = 0.0;
        double boundary_cuboid_max = 150.0;

        // creation of the mesh for ABM
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0u, false, boundary_cuboid_max-0, boundary_cuboid_max/2));
        nodes.push_back(new Node<2>(1u, false, boundary_cuboid_max-1, boundary_cuboid_max/2));
        nodes.push_back(new Node<2>(2u, false, boundary_cuboid_max-2, boundary_cuboid_max/2));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 50.0); // estimation of box domain 

        ///////////
        // CELLS // 
        ///////////

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide, fully constrained 
        MAKE_PTR(TransitCellProliferativeType, p_transit_type); // vessel segment ? 

        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);

        cells[0]->GetCellData()->SetItem("Connectivity", 0);
        cells[1]->GetCellData()->SetItem("Connectivity", 0);
        cells[2]->GetCellData()->SetItem("Connectivity", 1);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        cell_population.Update(); // addition of this line compared to the sprouting test

        // Writers 
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<ConsecutiveBranchesWriter>();
        cell_population.AddCellWriter<BranchNumberWriter>();
        cell_population.AddCellWriter<BirthTimeCellWriter>();
        cell_population.AddCellWriter<TortuosityWriter>();
        cell_population.AddCellWriter<AnastomosisWriter>();
        cell_population.AddCellWriter<ConnectivityWriter>();

        // fully constrain the first cell using the boundary condition 
        unsigned node_index_tip_cell = cell_population.GetLocationIndexUsingCell(0);
        std::vector<unsigned> pinned_node_indices;
        pinned_node_indices.push_back(node_index_tip_cell);
        MAKE_PTR_ARGS(PinnedCellsBoundaryCondition<2>, p_boundary_condition, (&cell_population, pinned_node_indices));

        //////////////
        // SETTINGS // 
        //////////////

        // set seed 
        RandomNumberGenerator::Instance()->Reseed(input_val_seed);

        // Set up simulation time for file output
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory_vegfconstant);
        simulator.SetSamplingTimestepMultiple(120*13);
        simulator.SetEndTime(input_val_time);
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

        /////////////////
        // SIMULATION // 
        ////////////////

        // 1) SOLVE PDE 

        // No PDE in the case of constant VEGF in the ECM

        // 2) UPDATING CELL POSITION 

        // Random force (all cells)
        MAKE_PTR_ARGS(RandomForce<2>, p_random_force, (input_val_sigma));
        simulator.AddForce(p_random_force);

        // Chemotactic force (tip cells only) 
        MAKE_PTR_ARGS(ChemoForceWithConstantVegf<2>, p_chemo_force, (input_val_chi, 1e-4)); 
        simulator.AddForce(p_chemo_force);

        //Persistence force (tip cells only)
        MAKE_PTR_ARGS(PersistenceForce<2>, p_persistence_force, (input_val_omegap));
        simulator.AddForce(p_persistence_force);

        // Mechanical force (all cells)
        MAKE_PTR(LinearMechanicalForceModified<2>, p_mechanical_force);
        p_mechanical_force->SetMeinekeSpringStiffness(15.0*input_val_S);
        cell_population.SetMeinekeDivisionSeparation(1.0);
        p_mechanical_force->SetMeinekeDivisionRestingSpringLength(1.0);
        p_mechanical_force->SetMeinekeSpringGrowthDuration(1.0);
        p_mechanical_force->SetCutOffLength(1.5);
        simulator.AddForce(p_mechanical_force);

        // Angular force (vessel segment only)
        MAKE_PTR_ARGS(AngularForce<2>, p_angular_force, (-input_val_omegaa)); 
        simulator.AddForce(p_angular_force);

        // 3) DIVISION OF CELLS 

        // Set the division rule for our population to be the random direction division rule
        MAKE_PTR_ARGS(SproutingRuleWithConstantVegf<2>, p_division_rule_to_set, (input_val_maxsproutingrate, input_val_thresholdlength, input_val_vegf_constantbackground, input_val_cmax, input_val_cmin, input_val_pmax, input_val_pmin));
        
        // Set the division rule for our population to be the new division rule implemented earlier 
        cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);

        // we set for each new daughter cell in the population if it is a tip cell or a vessel segment by using the function DaughterTypeofCell
        MAKE_PTR_ARGS(DaughterCellModifier<2>, p_daughtercell_modifier, (input_val_anastomosislength, input_val_thresholdlength));
        simulator.AddSimulationModifier(p_daughtercell_modifier);

        MAKE_PTR_ARGS(DirectionalPersistenceCellModifier<2>, p_persistenceforce_modifier, ());
        simulator.AddSimulationModifier(p_persistenceforce_modifier);

        cell_population.Update();

        simulator.Solve();

        // Output run time data
        CellBasedEventHandler::Headings();
        CellBasedEventHandler::Report();

        SimulationTime::Destroy();
    }


    void NoTestAngiogenesisModelWithVegfConcentrationConstantIn3D() 
    {
        ///////////
        // INPUT // 
        ///////////

        // Parameters input 
        CommandLineArguments* command_line = CommandLineArguments::Instance();

        // parameters for vegf concentration
        double input_val_vegf_constantbackground = command_line->GetDoubleCorrespondingToOption("-vegfconstantbackground"); //0.1

        // parameters for forces 
        double input_val_S = command_line->GetDoubleCorrespondingToOption("-S"); //1E-5
        double input_val_sigma = command_line->GetDoubleCorrespondingToOption("-sigma"); //1E-1
        double input_val_chi = command_line->GetDoubleCorrespondingToOption("-chi"); //1E-4
        double input_val_omegap = command_line->GetDoubleCorrespondingToOption("-omegap"); //1E-1
        double input_val_omegaa = command_line->GetDoubleCorrespondingToOption("-omegaa"); //1E-5

        // parameters for Psprout 
        double input_val_maxsproutingrate = command_line->GetDoubleCorrespondingToOption("-maxsproutingrate"); 
        double input_val_cmax = command_line->GetDoubleCorrespondingToOption("-cmax");
        double input_val_cmin = command_line->GetDoubleCorrespondingToOption("-cmin");
        double input_val_pmax = command_line->GetDoubleCorrespondingToOption("-pmax");
        double input_val_pmin = command_line->GetDoubleCorrespondingToOption("-pmin"); 

        // parameters for anastomosis 
        double input_val_anastomosislength = command_line->GetDoubleCorrespondingToOption("-anastomosislength");
        double input_val_thresholdlength = command_line->GetDoubleCorrespondingToOption("-thresholdlength");;

        // general parameters (time, random seed, output directory)
        double input_val_time = command_line->GetDoubleCorrespondingToOption("-time"); 
        double input_val_seed = command_line->GetIntCorrespondingToOption("-seed"); 
        std::string output_directory_constantvegf = command_line->GetStringCorrespondingToOption("-output_directory_vegfconstant");

        ////////////
        // MESHES // 
        ////////////

        // mesh for pdes 
        //double boundary_cuboid_min = 0.0;
        double boundary_cuboid_max = 150.0;

        // creation of the mesh for ABM
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u, false, boundary_cuboid_max-0, boundary_cuboid_max/2, boundary_cuboid_max/2));
        nodes.push_back(new Node<3>(1u, false, boundary_cuboid_max-1, boundary_cuboid_max/2, boundary_cuboid_max/2));
        nodes.push_back(new Node<3>(2u, false, boundary_cuboid_max-2, boundary_cuboid_max/2, boundary_cuboid_max/2));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 50.0); // estimation of box domain: go back to 1.5 or at least smaller than 50, try different numbers 

        ///////////
        // CELLS // 
        ///////////

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

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

        cells[0]->GetCellData()->SetItem("Connectivity", 0);
        cells[1]->GetCellData()->SetItem("Connectivity", 0);
        cells[2]->GetCellData()->SetItem("Connectivity", 1);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        cell_population.Update(); // addition of this line compared to the sprouting test

        // Writers 
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<ConsecutiveBranchesWriter>();
        cell_population.AddCellWriter<BranchNumberWriter>();
        cell_population.AddCellWriter<BirthTimeCellWriter>();
        cell_population.AddCellWriter<TortuosityWriter>();
        cell_population.AddCellWriter<AnastomosisWriter>();
        cell_population.AddCellWriter<ConnectivityWriter>();

        // fully constrain the first cell using the boundary condition 
        unsigned node_index_tip_cell = cell_population.GetLocationIndexUsingCell(0);
        std::vector<unsigned> pinned_node_indices;
        pinned_node_indices.push_back(node_index_tip_cell);
        MAKE_PTR_ARGS(PinnedCellsBoundaryCondition<3>, p_boundary_condition, (&cell_population, pinned_node_indices));

        //////////////
        // SETTINGS // 
        //////////////

        // set seed 
        RandomNumberGenerator::Instance()->Reseed(input_val_seed);

        // Set up simulation time for file output
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory_constantvegf);
        simulator.SetSamplingTimestepMultiple(120*13); // output every 100h ; cell cycle ~13 h 
        simulator.SetEndTime(input_val_time); // 2000
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

        /////////////////
        // SIMULATION // 
        ////////////////

        // 1) SOLVE PDE 

        // No PDE in the case of constant VEGF in the ECM 

        // 2) UPDATING CELL POSITION 

        // Random force (all cells)
        MAKE_PTR_ARGS(RandomForce<3>, p_random_force, (input_val_sigma));
        simulator.AddForce(p_random_force);

        // Chemotactic force (tip cells only) 
        MAKE_PTR_ARGS(ChemoForceWithConstantVegf<3>, p_chemo_force, (input_val_chi, 1e-4));
        simulator.AddForce(p_chemo_force);

        //Persistence force (tip cells only)
        MAKE_PTR_ARGS(PersistenceForce<3>, p_persistence_force, (input_val_omegap));
        simulator.AddForce(p_persistence_force);

        // Mechanical force (all cells)
        MAKE_PTR(LinearMechanicalForceModified<3>, p_mechanical_force);
        p_mechanical_force->SetMeinekeSpringStiffness(15.0*input_val_S);
        cell_population.SetMeinekeDivisionSeparation(1.0);
        p_mechanical_force->SetMeinekeDivisionRestingSpringLength(1.0);
        p_mechanical_force->SetMeinekeSpringGrowthDuration(1.0);
        p_mechanical_force->SetCutOffLength(1.5);
        simulator.AddForce(p_mechanical_force);

        // Angular force (vessel segment only)
        MAKE_PTR_ARGS(AngularForce<3>, p_angular_force, (-input_val_omegaa)); 
        simulator.AddForce(p_angular_force);

        // 3) DIVISION OF CELLS 

        // Set the division rule for our population to be the random direction division rule
        MAKE_PTR_ARGS(SproutingRuleWithConstantVegf<3>, p_division_rule_to_set, (input_val_maxsproutingrate, input_val_thresholdlength, input_val_vegf_constantbackground, input_val_cmax, input_val_cmin, input_val_pmax, input_val_pmin));

        // Set the division rule for our population to be the new division rule implemented earlier 
        cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);

        // we set for each new daughter cell in the population if it is a tip cell or a vessel segment by using the function DaughterTypeofCell
        MAKE_PTR_ARGS(DaughterCellModifier<3>, p_daughtercell_modifier, (input_val_anastomosislength, input_val_thresholdlength));
        simulator.AddSimulationModifier(p_daughtercell_modifier);

        MAKE_PTR_ARGS(DirectionalPersistenceCellModifier<3>, p_persistenceforce_modifier, ());
        simulator.AddSimulationModifier(p_persistenceforce_modifier);

        cell_population.Update();

        simulator.Solve();

        // Output run time data
        CellBasedEventHandler::Headings();
        CellBasedEventHandler::Report();

        SimulationTime::Destroy();
    }

};

#endif /*TESTANGIOGENESISMODELWITHVEGFCONCENTRATIONCONSTANT_HPP_*/