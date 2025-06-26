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

    // we test each forces independently first, without division rule

    void TestRandomForce()
    {
        // creation of the mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0u, false, boundary_cuboid_max-20, boundary_cuboid_max/2));
        nodes.push_back(new Node<2>(1u, false, boundary_cuboid_max-21, boundary_cuboid_max/2));
        nodes.push_back(new Node<2>(2u, false, boundary_cuboid_max-22, boundary_cuboid_max/2));

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

        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        // we create the random force object 
        RandomForce<2> random_force(0.1);

        random_force.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 0.1, 1e-4); // x axis
            TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 0.1, 1e-4); // y axis 
        }
    }

    void NoTestChemotacticForce()
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

    void NoTestPersistenceForce()
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

    void NoTestAngularForce() throw(Exception)
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

    void NoTestSproutingRule() 
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
        MAKE_PTR_ARGS(DaughterCellModifier<2>, p_modifier, ());
        simulator.AddSimulationModifier(p_modifier);

        // run simulation 
        simulator.Solve();

        SimulationTime::Destroy();
    }

};

#endif /*TESTFORCESMODEL_HPP_*/