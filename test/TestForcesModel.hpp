#ifndef TESTFORCESMODEL_HPP_
#define TESTFORCESMODEL_HPP_

#include <cxxtest/TestSuite.h>

#include "CellPropertyRegistry.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "OutputFileHandler.hpp"

#include "CheckpointArchiveTypes.hpp"

#include "GeneralisedLinearSpringForce.hpp"
#include "DifferentialAdhesionGeneralisedLinearSpringForce.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "CellLabel.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "OffLatticeSimulation.hpp"
#include "CellData.hpp"
#include "CellId.hpp"
#include "UniformCellCycleModel.hpp"
#include "AbstractCellPopulationBoundaryCondition.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "AbstractCellMutationState.hpp"
#include "Cell.hpp"
#include "RandomDirectionCentreBasedDivisionRule.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellBasedEventHandler.hpp"

#include "PetscSetupAndFinalize.hpp"
#include "Debug.hpp"

#include "ChemoForce.hpp"
#include "RandomForce.hpp"
#include "PersistenceForce.hpp"
#include "AngularForce.hpp"
#include "LinearMechanicalForceModified.hpp"
#include "DirectionalPersistenceCellModifier.hpp"

#include "SproutingRule.hpp"
#include "DaughterCellModifier.hpp"
#include "PinnedCellsBoundaryCondition.hpp"

#include "BranchingCellMutationState.hpp"
#include "TipCellMutationState.hpp"
#include "VesselCellMutationState.hpp"

class TestForcesModel : public AbstractCellBasedTestSuite
{

public:

    // we test each forces independently first, without division rule and with our cell cycle model

    // RANDOM FORCE TEST //
    void TestRandomForce()
        {
            EXIT_IF_PARALLEL; // Honeycomb mesh not made to be run in parallel

            // creation of the mesh
            HoneycombMeshGenerator generator(3,1);
            MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal'sw model)

            // creation of the cells 
            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

            // creation of a population of cells 
            NodeBasedCellPopulation<2> cell_population(mesh, cells);

            // creation of the simulation 
            SimulationTime::Instance() -> SetStartTime(0.0);
            OffLatticeSimulation<2> simulator(cell_population);
            simulator.SetOutputDirectory("TestRandomForce");
            simulator.SetSamplingTimestepMultiple(12);
            simulator.SetEndTime(5.0);

            // we apply the random force law 
            MAKE_PTR_ARGS(RandomForce<2>, p_random_force, (0.4));
            simulator.AddForce(p_random_force);

            // run simulation 
            simulator.Solve();

            SimulationTime::Destroy();
        }

    // CHEMOTACTIC FORCE TEST //
    void NTestChemoForce()
        {
            EXIT_IF_PARALLEL; // Honeycomb mesh not made to be run in parallel

            // creation of the mesh
            HoneycombMeshGenerator generator(3,1);
            MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal'sw model)

            // creation of the cells 
            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

            // creation of a population of cells 
            NodeBasedCellPopulation<2> cell_population(mesh, cells);

            // creation of the simulation 
            SimulationTime::Instance() -> SetStartTime(0.0);
            OffLatticeSimulation<2> simulator(cell_population);
            simulator.SetOutputDirectory("TestChemoForce");
            simulator.SetSamplingTimestepMultiple(12);
            simulator.SetEndTime(5.0);

            // we apply the random force law 
            MAKE_PTR_ARGS(ChemoForce<2>, p_chemo_force, (0.1, 5.56));
            simulator.AddForce(p_chemo_force);

            // run simulation 
            simulator.Solve();

            SimulationTime::Destroy();
        }

        // PERSISTENCE FORCE // 
        void NTestPersistenceForce()
        {
            EXIT_IF_PARALLEL; // Honeycomb mesh not made to be run in parallel

            // creation of the mesh
            HoneycombMeshGenerator generator(3,1);
            MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal'sw model)

            // creation of the cells 
            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

            // creation of a population of cells 
            NodeBasedCellPopulation<2> cell_population(mesh, cells);
            cell_population.Update();

            // creation of the simulation 
            SimulationTime::Instance() -> SetStartTime(0.0);
            OffLatticeSimulation<2> simulator(cell_population);
            simulator.SetOutputDirectory("TestPersistenceForce");
            simulator.SetSamplingTimestepMultiple(12);
            simulator.SetEndTime(5.0);

            for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
            {
                cell_population.GetNode(i)->ClearAppliedForce();
            }

            // MAKE_PTR_ARGS(GeneralisedLinearSpringForce<2>, p_linear_force, ());
            // simulator.AddForce(p_linear_force);
            MAKE_PTR_ARGS(ChemoForce<2>, p_chemo_force, (0.1, 5.56));
            simulator.AddForce(p_chemo_force);

            unsigned node_index = 0;
            // Set up cell data on the cell population
            for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
            cell_iter != cell_population.End();
            ++cell_iter) 
            {
                CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
                c_vector<double, 2> new_r_cellmovement = cell_population.GetLocationOfCellCentre(p_cell);
                cell_iter->GetCellData()-> SetItem("old_x_coordinate", new_r_cellmovement(0));
                cell_iter->GetCellData()-> SetItem("old_y_coordinate", new_r_cellmovement(1));
                ++node_index;
            }

            // we apply the persistence force law 
            MAKE_PTR_ARGS(PersistenceForce<2>, p_persistence_force, (0.4));
            simulator.AddForce(p_persistence_force);

            // run simulation 
            simulator.Solve();

            SimulationTime::Destroy();
        }

        // ANGULAR FORCE // 
        void NTestAngularForce() throw(Exception)
        {
            EXIT_IF_PARALLEL; // Honeycomb mesh not made to be run in parallel

            // creation of the mesh
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0u, false, 1.0, 0.0));
            nodes.push_back(new Node<2>(1u, false, 1.0, 1.5));
            nodes.push_back(new Node<2>(2u, false, 1.5, 0.6));
            nodes.push_back(new Node<2>(3u, false, 0.6, 1.0));

            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal'sw model)

            // creation of the cells 
            std::vector<CellPtr> cells;
            MAKE_PTR(StemCellProliferativeType, p_stem_type);
            MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);
            cells[2]->SetCellProliferativeType(p_stem_type);
            cells[3]->SetCellProliferativeType(p_stem_type);

            // creation of a population of cells 
            NodeBasedCellPopulation<2> cell_population(mesh, cells);
            cell_population.Update();

            // initialisation of the forces 
            for(unsigned i=0; i<cell_population.GetNumNodes(); i++)
            {
                cell_population.GetNode(i)->ClearAppliedForce();
            }

            // creation of the simulation 
            SimulationTime::Instance() -> SetStartTime(0.0);
            OffLatticeSimulation<2> simulator(cell_population);
            simulator.SetOutputDirectory("TestAngularForce");
            simulator.SetSamplingTimestepMultiple(12);
            simulator.SetEndTime(5.0);

            // we apply the angular force law 
            MAKE_PTR_ARGS(AngularForce<2>, p_angular_force, (5.56));
            simulator.AddForce(p_angular_force);

            // run simulation 
            simulator.Solve();

            SimulationTime::Destroy();
        }

        // we test now the sprouting rule 
       void NTestSproutingRule() 
       {
           EXIT_IF_PARALLEL; // Honeycomb mesh not made to be run in parallel

           // creation of the mesh
           HoneycombMeshGenerator generator(3,1);
           MutableMesh<2,2>* p_generating_mesh = generator.GetMesh(); //Mutable and not Tetrahedral Mesh compared to the Sprouting test
           NodesOnlyMesh<2> mesh;
           mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal'sw model)

           // creation of the cells 
           std::vector<CellPtr> cells;
           MAKE_PTR(StemCellProliferativeType, p_stem_type);
           CellsGenerator<UniformCellCycleModel, 2> cells_generator;
           cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_stem_type);

           // creation of a population of cells 
           NodeBasedCellPopulation<2> cell_population(mesh, cells);
           cell_population.Update(); // addition of this line compared to the sprouting test

           // we copy this cell population to obtain all the previous location index from the cell population before the division 
           unsigned old_number_nodes = cell_population.GetNumNodes();

           // Set the division rule for our population to be the random direction division rule
           typedef SproutingRule<2,2> SproutingRule;
           MAKE_PTR_ARGS(SproutingRule, p_division_rule_to_set, ());

           // initialisation of the forces : not in the sprouting test
           for(unsigned i=0; i<cell_population.GetNumNodes(); i++)
           {
               cell_population.GetNode(i)->ClearAppliedForce();
           }

           OffLatticeSimulation<2> simulator(cell_population);
           simulator.SetOutputDirectory("TestSproutingRule");
           simulator.SetSamplingTimestepMultiple(12);
           simulator.SetEndTime(20.0);

           // Set the division rule for our population to be the new division rule implemented earlier 
           //boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule_to_set(new SproutingRule()); // do not function 
           cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);
           //SproutingRule<2,2>* p_sprouting_division_rule; // do not function 

            // we set for each new daughter cell in the population if it is a tip cell or a vessel segment by using the function DaughterTypeofCell
            MAKE_PTR_ARGS(DaughterCellModifier<2>, p_modifier, (old_number_nodes));
            simulator.AddSimulationModifier(p_modifier);

           // run simulation 
           simulator.Solve();

           SimulationTime::Destroy();
       }

    // we test the sprouting division rule with the cells property : tip cell - vessel segment - tip cell 
       void NTestSproutingRuleAndCellCycle() 
       {
            EXIT_IF_PARALLEL; // Honeycomb mesh not made to be run in parallel

            HoneycombMeshGenerator generator(3,1);
            std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
            MutableMesh<2,2>* p_generating_mesh = generator.GetMesh(); //Mutable and not Tetrahedral Mesh compared to the Sprouting test
            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal'sw model)

            // creation of the cells 
            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // tip cells
            MAKE_PTR(StemCellProliferativeType, p_stem_type); // vessel segment
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            cells_generator.GenerateGivenLocationIndices(cells, location_indices, p_differentiated_type);
            cells[1]->SetCellProliferativeType(p_stem_type);
           
           // creation of a population of cells 
           NodeBasedCellPopulation<2> cell_population(mesh, cells);
           cell_population.Update(); // addition of this line compared to the sprouting test

           // we copy this cell population to obtain all the previous location index from the cell population before the division 
           unsigned old_number_nodes = cell_population.GetNumNodes();

           // Set the division rule for our population to be the random direction division rule
           typedef SproutingRule<2,2> SproutingRule;
           MAKE_PTR_ARGS(SproutingRule, p_division_rule_to_set, ());

           // initialisation of the forces : not in the sprouting test
           for(unsigned i=0; i<cell_population.GetNumNodes(); i++)
           {
               cell_population.GetNode(i)->ClearAppliedForce();
           }

           OffLatticeSimulation<2> simulator(cell_population);
           simulator.SetOutputDirectory("TestSproutingRuleAndCellCycle");
           simulator.SetSamplingTimestepMultiple(12);
           simulator.SetEndTime(20.0);

           // Set the division rule for our population to be the new division rule implemented earlier 
           //boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule_to_set(new SproutingRule()); // do not function 
           cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);
           //SproutingRule<2,2>* p_sprouting_division_rule; // do not function 

           // we set for each new daughter cell in the population if it is a tip cell or a vessel segment by using the function DaughterTypeofCell
            MAKE_PTR_ARGS(DaughterCellModifier<2>, p_modifier, (old_number_nodes));
            simulator.AddSimulationModifier(p_modifier);

           // run simulation 
           simulator.Solve();

           cout << "Size of new cell population = " << cell_population.GetNumNodes() << endl;

           SimulationTime::Destroy();
       }

    void NTestMechanicalForceAndAngularForceForVesselSegments() 
       {
            EXIT_IF_PARALLEL; // Honeycomb mesh not made to be run in parallel
            HoneycombMeshGenerator generator(3,1);

            // creation of the mesh
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0u, false, 0.0, 0.0));
            nodes.push_back(new Node<2>(1u, false, 0.8, 0.0));
            nodes.push_back(new Node<2>(2u, false, 1.55, 0.0));
            std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal'sw model)

            // creation of the cells 
            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // tip cell 
            MAKE_PTR(StemCellProliferativeType, p_stem_type); // vessel segment
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            cells_generator.GenerateGivenLocationIndices(cells, location_indices, p_differentiated_type);
            cells[1]->SetCellProliferativeType(p_stem_type);
           
           // creation of a population of cells 
           NodeBasedCellPopulation<2> cell_population(mesh, cells);
           cell_population.Update(); // addition of this line compared to the sprouting test

           // we copy this cell population to obtain all the previous location index from the cell population before the division 
           unsigned OldNumNodes = cell_population.GetNumNodes();

           // Set up cell data on the cell population : initialisation for the division 
            for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
            cell_iter != cell_population.End();
            ++cell_iter)
            {
                unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);
                CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
                cell_iter->GetCellData()->SetItem("node_index", node_index);
                cell_iter->GetCellData()->SetItem("daughter_type", 0.0);
            }

           // initialisation of the forces : not in the sprouting test
           for(unsigned i=0; i<cell_population.GetNumNodes(); i++)
           {
               cell_population.GetNode(i)->ClearAppliedForce();
           }

           OffLatticeSimulation<2> simulator(cell_population);
           simulator.SetOutputDirectory("TestMechanicalForceAndAngularForceForVesselSegments");
           simulator.SetSamplingTimestepMultiple(12);
           simulator.SetEndTime(50.0);

            /////////////////
            // SIMULATION // 
            ////////////////

            // 1) UPDATING CELL POSITION 

            // Mechanical force (all cells)
            MAKE_PTR(GeneralisedLinearSpringForce<2>, p_mechanical_force);
            //p_mechanical_force->SetMeinekeSpringStiffness(5.56);
            simulator.AddForce(p_mechanical_force);

            // Angular force (vessel segment only)
            MAKE_PTR_ARGS(AngularForce<2>, p_angular_force, (5.56E-3));
            simulator.AddForce(p_angular_force);

            // 2) DIVISION OF CELLS 

            // Set the division rule for our population to be the random direction division rule
           typedef SproutingRule<2,2> SproutingRule;

           // Set the division rule for our population to be the new division rule implemented earlier 
           MAKE_PTR_ARGS(SproutingRule, p_division_rule_to_set, (0.5));
           cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);

            // we set for each new daughter cell in the population if it is a tip cell or a vessel segment by using the function DaughterTypeofCell
            MAKE_PTR_ARGS(DaughterCellModifier<2>, p_modifier, (OldNumNodes));
            simulator.AddSimulationModifier(p_modifier);

            cell_population.Update();

            simulator.Solve();

           SimulationTime::Destroy();
       }

};

#endif /*TESTFORCESMODEL_HPP_*/