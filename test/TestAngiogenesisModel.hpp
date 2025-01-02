#ifndef TESTANGIOGENESISMODEL_HPP_
#define TESTANGIOGENESISMODEL_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>
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

// Cell properties
#include "CellPropertyRegistry.hpp"
#include "CellLabel.hpp"
#include "CellData.hpp"
#include "CellId.hpp"
#include "Cell.hpp"
#include "CellBasedEventHandler.hpp"

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

// Forces 
#include "GeneralisedLinearSpringForce.hpp"
#include "DifferentialAdhesionGeneralisedLinearSpringForce.hpp"
#include "ChemoForce.hpp"
#include "RandomForce.hpp"
#include "PersistenceForce.hpp"
#include "AngularForce.hpp"
#include "LinearMechanicalForceModified.hpp"
#include "DirectionalPersistenceCellModifier.hpp"

// Writers
#include "CellMutationStatesCountWriter.hpp"
#include "ConsecutiveBranchesWriter.hpp"
#include "BranchesNumberWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "BirthTimeCellWriter.hpp"
#include "TortuosityWriter.hpp"


// PDE solvers
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "SimpleNonlinearEllipticSolver.hpp"
#include "FunctionalBoundaryCondition.hpp"
#include "SimpleNewtonNonlinearSolver.hpp"
#include "VegfEquationPde.hpp"
#include "TrianglesMeshReader.hpp"


class TestAngiogenesisModel : public AbstractCellBasedTestSuite
{

public:

    // we test the different forces with the cell cycle and the new division rule 
    // we need to check for each cell if it is a tip cell or a vessel segment 
    void NoTestSproutingRuleAndCellCycleAndForcesForAllCellsIN2D() 
       {
            // to change the values of the test directly on the command line 
            CommandLineArguments* command_line = CommandLineArguments::Instance();
            double input_val_sigma = command_line->GetDoubleCorrespondingToOption("-sigma");
            double input_val_chi = command_line->GetDoubleCorrespondingToOption("-chi");
            double input_val_omegap = command_line->GetDoubleCorrespondingToOption("-omegap");
            double input_val_omegaa = command_line->GetDoubleCorrespondingToOption("-omegaa");
            double input_val_Psprout = command_line->GetDoubleCorrespondingToOption("-Psprout");
            double input_val_time = command_line->GetDoubleCorrespondingToOption("-time");
            double input_val_seed = command_line->GetIntCorrespondingToOption("-seed");

            std::string output_directory = command_line->GetStringCorrespondingToOption("-output_directory");

            // set seed 
            RandomNumberGenerator::Instance()->Reseed(input_val_seed);

            // creation of the mesh
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0u, false, 0.0, 0.0));
            nodes.push_back(new Node<2>(1u, false, 1.0, 0.0));
            nodes.push_back(new Node<2>(2u, false, 2.0, 0.0));

            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal'sw model)

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
            cell_population.AddCellWriter<ConsecutiveBranchesWriter>();
            cell_population.AddCellWriter<BranchesNumberWriter>();
            cell_population.AddCellWriter<BirthTimeCellWriter>();
            cell_population.AddCellWriter<TortuosityWriter>();

            // we copy this cell population to obtain all the previous location index from the cell population before the division 
            unsigned OldNumNodes = cell_population.GetNumNodes();

            unsigned node_index_tip_cell = cell_population.GetLocationIndexUsingCell(0);

            // fully constrain the first cell using the boundary condition 
            std::vector<unsigned> pinned_node_indices;
            pinned_node_indices.push_back(node_index_tip_cell);

            typedef PinnedCellsBoundaryCondition<2,2> PinnedCellsBoundaryCondition;
            MAKE_PTR_ARGS(PinnedCellsBoundaryCondition, p_boundary_condition, (&cell_population, pinned_node_indices));

            OffLatticeSimulation<2> simulator(cell_population);
            simulator.SetOutputDirectory(output_directory);
            simulator.SetSamplingTimestepMultiple(120);
            simulator.SetEndTime(input_val_time);
            simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

            /////////////////
            // SIMULATION // 
            ////////////////

            // 1) UPDATING CELL POSITION 

            // Random force (all cells)
            MAKE_PTR_ARGS(RandomForce<2>, p_random_force, (input_val_sigma));
            simulator.AddForce(p_random_force);

            // Chemotactic force (tip cells only) 
            MAKE_PTR_ARGS(ChemoForce<2>, p_chemo_force, (-input_val_chi, 1E-2, 0.0));
            simulator.AddForce(p_chemo_force);

            //Persistence force (tip cells only)
            MAKE_PTR_ARGS(PersistenceForce<2>, p_persistence_force, (input_val_omegap));
            simulator.AddForce(p_persistence_force);

            // Mechanical force (all cells)
            MAKE_PTR(LinearMechanicalForceModified<2>, p_mechanical_force);
            p_mechanical_force->SetMeinekeSpringStiffness(15.0);
            cell_population.SetMeinekeDivisionSeparation(1.0);
            p_mechanical_force->SetMeinekeDivisionRestingSpringLength(1.0);
            p_mechanical_force->SetMeinekeSpringGrowthDuration(1.0);
            p_mechanical_force->SetCutOffLength(1.5);
            simulator.AddForce(p_mechanical_force);

            // Angular force (vessel segment only)
            MAKE_PTR_ARGS(AngularForce<2>, p_angular_force, (-input_val_omegaa));
            simulator.AddForce(p_angular_force);

            // 2) DIVISION OF CELLS 

            // Set the division rule for our population to be the random direction division rule
            typedef SproutingRule<2,2> SproutingRule;
            MAKE_PTR_ARGS(SproutingRule, p_division_rule_to_set, (input_val_Psprout));

            // Set the division rule for our population to be the new division rule implemented earlier 
            cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);

            // we set for each new daughter cell in the population if it is a tip cell or a vessel segment by using the function DaughterTypeofCell
            MAKE_PTR_ARGS(DaughterCellModifier<2>, p_modifier, (OldNumNodes));
            simulator.AddSimulationModifier(p_modifier);

            MAKE_PTR_ARGS(DirectionalPersistenceCellModifier<2>, p_modifier_2, ());
            simulator.AddSimulationModifier(p_modifier_2);

            cell_population.Update();

            simulator.Solve();

            // Output run time data
            CellBasedEventHandler::Headings();
            CellBasedEventHandler::Report();

           SimulationTime::Destroy();
       }

       void TestSproutingRuleAndCellCycleAndForcesForAllCellsIN3D() 
       {
            // to change the values of the test directly on the command line 
            CommandLineArguments* command_line = CommandLineArguments::Instance();
            double input_val_sigma = command_line->GetDoubleCorrespondingToOption("-sigma");
            double input_val_chi = command_line->GetDoubleCorrespondingToOption("-chi");
            double input_val_omegap = command_line->GetDoubleCorrespondingToOption("-omegap");
            double input_val_omegaa = command_line->GetDoubleCorrespondingToOption("-omegaa");
            double input_val_Psprout = command_line->GetDoubleCorrespondingToOption("-Psprout");
            double input_val_time = command_line->GetDoubleCorrespondingToOption("-time");
            double input_val_seed = command_line->GetIntCorrespondingToOption("-seed");

            std::string output_directory = command_line->GetStringCorrespondingToOption("-output_directory");

            // set seed 
            RandomNumberGenerator::Instance()->Reseed(input_val_seed);

            // creation of the mesh
            std::vector<Node<3>*> nodes;
            nodes.push_back(new Node<3>(0u, false, 0.0, 0.0, 0.0));
            nodes.push_back(new Node<3>(1u, false, 1.0, 0.0, 0.0));
            nodes.push_back(new Node<3>(2u, false, 2.0, 0.0, 0.0));

            NodesOnlyMesh<3> mesh;
            mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

            // initialisation of PDEs and their solvers

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
            cell_population.AddCellWriter<ConsecutiveBranchesWriter>();
            cell_population.AddCellWriter<BranchesNumberWriter>();
            cell_population.AddCellWriter<BirthTimeCellWriter>();
            cell_population.AddCellWriter<TortuosityWriter>();

            // we copy this cell population to obtain all the previous location index from the cell population before the division 
            unsigned OldNumNodes = cell_population.GetNumNodes();

            unsigned node_index_tip_cell = cell_population.GetLocationIndexUsingCell(0);

            // fully constrain the first cell using the boundary condition 
            std::vector<unsigned> pinned_node_indices;
            pinned_node_indices.push_back(node_index_tip_cell);

            typedef PinnedCellsBoundaryCondition<3,3> PinnedCellsBoundaryCondition;
            MAKE_PTR_ARGS(PinnedCellsBoundaryCondition, p_boundary_condition, (&cell_population, pinned_node_indices));

            OffLatticeSimulation<3> simulator(cell_population);
            simulator.SetOutputDirectory(output_directory);
            simulator.SetSamplingTimestepMultiple(120);
            simulator.SetEndTime(input_val_time);
            simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

            /////////////////
            // SIMULATION // 
            ////////////////

            // 1) UPDATING CELL POSITION 

            // Random force (all cells)
            MAKE_PTR_ARGS(RandomForce<3>, p_random_force, (input_val_sigma));
            simulator.AddForce(p_random_force);

            // Chemotactic force (tip cells only) 
            MAKE_PTR_ARGS(ChemoForce<3>, p_chemo_force, (-input_val_chi, 1.0E-2, 0.0));
            simulator.AddForce(p_chemo_force);

            //Persistence force (tip cells only)
            MAKE_PTR_ARGS(PersistenceForce<3>, p_persistence_force, (input_val_omegap));
            simulator.AddForce(p_persistence_force);

            // Mechanical force (all cells)
            MAKE_PTR(LinearMechanicalForceModified<3>, p_mechanical_force);
            p_mechanical_force->SetMeinekeSpringStiffness(15.0);
            cell_population.SetMeinekeDivisionSeparation(1.0);
            p_mechanical_force->SetMeinekeDivisionRestingSpringLength(1.0);
            p_mechanical_force->SetMeinekeSpringGrowthDuration(1.0);
            p_mechanical_force->SetCutOffLength(1.5);
            simulator.AddForce(p_mechanical_force);

            // Angular force (vessel segment only)
            MAKE_PTR_ARGS(AngularForce<3>, p_angular_force, (-input_val_omegaa));
            simulator.AddForce(p_angular_force);

            // 2) DIVISION OF CELLS 

            // Set the division rule for our population to be the random direction division rule
            typedef SproutingRule<3,3> SproutingRule;
            MAKE_PTR_ARGS(SproutingRule, p_division_rule_to_set, (input_val_Psprout));

            // Set the division rule for our population to be the new division rule implemented earlier 
            cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);

            // we set for each new daughter cell in the population if it is a tip cell or a vessel segment by using the function DaughterTypeofCell
            MAKE_PTR_ARGS(DaughterCellModifier<3>, p_modifier, (OldNumNodes));
            simulator.AddSimulationModifier(p_modifier);

            MAKE_PTR_ARGS(DirectionalPersistenceCellModifier<3>, p_modifier_2, ());
            simulator.AddSimulationModifier(p_modifier_2);

            cell_population.Update();

            simulator.Solve();

            // PDE solver (VEGF and MMP concentrations) 
            // TO WRITE

            // Output run time data
            CellBasedEventHandler::Headings();
            CellBasedEventHandler::Report();

            SimulationTime::Destroy();
       }

};

#endif /*TESTANGIOGENESISMODEL_HPP_*/