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
#include "ChemoForceWithConstantVegf.hpp"
#include "ChemoForceWithAnalyticalApproximationPde.hpp"
#include "ChemoForceWithPdes.hpp"
#include "RandomForce.hpp"
#include "PersistenceForce.hpp"
#include "AngularForce.hpp"
#include "LinearMechanicalForceModified.hpp"
#include "DirectionalPersistenceCellModifier.hpp"

#include "SproutingRule.hpp"
#include "SproutingRuleWithConstantVegf.hpp"
#include "SproutingRuleWithAnalyticalApproximationPde.hpp"
#include "SproutingRuleWithPdes.hpp"
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
        nodes.push_back(new Node<2>(0u, false, 10, 5));
        nodes.push_back(new Node<2>(1u, false, 9, 5));
        nodes.push_back(new Node<2>(2u, false, 8, 5));

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
            TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 0.1, 1e-2);
            TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 0.1, 1e-2);
            TS_ASSERT_EQUALS(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.1); // random norm
        }

        // test archiving class 
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "random_force.arch";
        {
            AbstractForce<2>* const p_force = new RandomForce<2>(0.1);
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_force;
            delete p_force;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractForce<2>* p_force;
            input_arch >> p_force;

            TS_ASSERT_EQUALS(dynamic_cast<RandomForce<2>*>(p_force)->GetRandomSensitivity(), 0.1);

            delete p_force;
        }
    }

    void TestMechanicalForce()
    {
        // creation of the mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0u, false, 10, 5));
        nodes.push_back(new Node<2>(1u, false, 9, 5));
        nodes.push_back(new Node<2>(2u, false, 8, 5));

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
        LinearMechanicalForceModified<2> mechanical_force;

        mechanical_force.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 1e-1, 1e-1); // random norm
        }

        // test archiving class 
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "random_force.arch";
        {
            AbstractForce<2>* const p_force = new LinearMechanicalForceModified<2>();
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_force;
            delete p_force;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractForce<2>* p_force;
            input_arch >> p_force;

            delete p_force;
        }
    }

    void TestChemotacticForce()
    {
        // creation of the mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0u, false, 10, 5));
        nodes.push_back(new Node<2>(1u, false, 9, 5));
        nodes.push_back(new Node<2>(2u, false, 8, 5));

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
        ChemoForce<2> chemo_force(1e-4, 0.1, -0.1, 0.0);

        chemo_force.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[0], -1e-5);
            TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 1e-5);
            TS_ASSERT_EQUALS(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 1e-5); // random norm
        }

        // test archiving class 
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "chemo_force.arch";
        {
            AbstractForce<2>* const p_force = new ChemoForce<2>(1e-4);
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_force;
            delete p_force;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractForce<2>* p_force;
            input_arch >> p_force;

            TS_ASSERT_EQUALS(dynamic_cast<ChemoForce<2>*>(p_force)->GetChemotacticSensitivity(), 1e-4);

            delete p_force;
        }
    }

    void TestChemotacticForceWithConstantVegf()
    {
        // creation of the mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0u, false, 10, 5));
        nodes.push_back(new Node<2>(1u, false, 9, 5));
        nodes.push_back(new Node<2>(2u, false, 8, 5));

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
        ChemoForceWithConstantVegf<2> chemo_force_withconstantvegf(1e-4, 0.1);

        chemo_force_withconstantvegf.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[0], -1e-5);
            TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 0.0);
            TS_ASSERT_EQUALS(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 1e-5); // random norm
        }

        // test archiving class 
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "chemo_force_withconstantvegf.arch";
        {
            AbstractForce<2>* const p_force = new ChemoForceWithConstantVegf<2>(1e-4);
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_force;
            delete p_force;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractForce<2>* p_force;
            input_arch >> p_force;

            TS_ASSERT_EQUALS(dynamic_cast<ChemoForceWithConstantVegf<2>*>(p_force)->GetChemotacticSensitivity(), 1e-4);

            delete p_force;
        }
    }

    void NoTestChemotacticForceWithPde() // to complete when study of PDE
    {
    }

    void TestChemotacticForceWithApproximationOfPde()
    {
        // creation of the mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0u, false, 10, 5));
        nodes.push_back(new Node<2>(1u, false, 9, 5));
        nodes.push_back(new Node<2>(2u, false, 8, 5));

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
        ChemoForceWithAnalyticalApproximationPde<2> chemo_force_withanalyticalapproximationpde(1e-4, 1e4, 1.0, 0.0, 0.0, 0.5, 0.1);

        chemo_force_withanalyticalapproximationpde.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            // choose one point and calculate the value at this point for the chemotactic force 
            TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[0], -3.7e-3, 1e-4); // -3.7e-3
            TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 0.0);
            TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 1e-5, 1e-6); 
        }

        // test archiving class 
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "chemo_force_withanalyticalapproximationpde.arch";
        {
            AbstractForce<2>* const p_force = new ChemoForceWithAnalyticalApproximationPde<2>(1e-4);
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_force;
            delete p_force;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractForce<2>* p_force;
            input_arch >> p_force;

            TS_ASSERT_EQUALS(dynamic_cast<ChemoForceWithAnalyticalApproximationPde<2>*>(p_force)->GetChemotacticSensitivity(), 1e-4);

            delete p_force;
        }
    }

    void TestPersistenceForce()
    {
        // creation of the mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0u, false, 10, 5));
        nodes.push_back(new Node<2>(1u, false, 9, 5));
        nodes.push_back(new Node<2>(2u, false, 8, 5));

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
        PersistenceForce<2> persistence_force(0.1);

        persistence_force.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 0.1, 1e-2);
            TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 0.1, 1e-2);
            TS_ASSERT_EQUALS(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.1); // random norm
        }

        // test archiving class 
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "persistence_force.arch";
        {
            AbstractForce<2>* const p_force = new PersistenceForce<2>(0.1);
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_force;
            delete p_force;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractForce<2>* p_force;
            input_arch >> p_force;

            TS_ASSERT_EQUALS(dynamic_cast<PersistenceForce<2>*>(p_force)->GetPersistenceCoefficient(), 0.1);

            delete p_force;
        }
    }

    void TestAngularForce()
    {
        // creation of the mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0u, false, 10, 5));
        nodes.push_back(new Node<2>(1u, false, 9, 5));
        nodes.push_back(new Node<2>(2u, false, 8, 5));

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
        AngularForce<2> angular_force(1e-5);

        angular_force.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            TS_ASSERT_EQUALS(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 1e-5); // random norm
        }

        // test archiving class 
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "angular_force.arch";
        {
            AbstractForce<2>* const p_force = new AngularForce<2>(1e-5);
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_force;
            delete p_force;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractForce<2>* p_force;
            input_arch >> p_force;

            TS_ASSERT_EQUALS(dynamic_cast<AngularForce<2>*>(p_force)->GetAngularPersistence(), 1e-5);

            delete p_force;
        }
    }

    void TestCellMutationStates()
    {
        // creation of the mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0u, false, 10, 5)); // pinned cell 
        nodes.push_back(new Node<2>(1u, false, 9, 5)); // stalk cell 
        nodes.push_back(new Node<2>(2u, false, 8, 5)); // branching cell 
        nodes.push_back(new Node<2>(3u, false, 7, 5.5)); // tip cell 1
        nodes.push_back(new Node<2>(4u, false, 7, 4.5)); // tip cell 2

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingCellMutationState, p_branching_state); 
        MAKE_PTR(VesselCellMutationState, p_vessel_state);
        MAKE_PTR(TipCellMutationState, p_tip_state);

        // test the pointers: colours and type 
        TS_ASSERT_EQUALS(p_branching_state->GetColour(), 2u);
        TS_ASSERT_EQUALS(p_branching_state->IsType<TipCellMutationState>(), false);
        TS_ASSERT_EQUALS(p_branching_state->IsType<VesselCellMutationState>(), false);
        TS_ASSERT_EQUALS(p_branching_state->IsType<BranchingCellMutationState>(), true);

        TS_ASSERT_EQUALS(p_vessel_state->GetColour(), 1u);
        TS_ASSERT_EQUALS(p_vessel_state->IsType<TipCellMutationState>(), false);
        TS_ASSERT_EQUALS(p_vessel_state->IsType<VesselCellMutationState>(), true);
        TS_ASSERT_EQUALS(p_vessel_state->IsType<BranchingCellMutationState>(), false);

        TS_ASSERT_EQUALS(p_tip_state->GetColour(), 0u);
        TS_ASSERT_EQUALS(p_tip_state->IsType<TipCellMutationState>(), true);
        TS_ASSERT_EQUALS(p_tip_state->IsType<VesselCellMutationState>(), false);
        TS_ASSERT_EQUALS(p_tip_state->IsType<BranchingCellMutationState>(), false);

        // test the number of cells and attribution of mutation state 
        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_branching_state);
        cells[3]->SetMutationState(p_tip_state);
        cells[4]->SetMutationState(p_tip_state);

        TS_ASSERT_EQUALS(p_branching_state->GetCellCount(), 1u);
        TS_ASSERT_EQUALS(p_vessel_state->GetCellCount(), 2u);
        TS_ASSERT_EQUALS(p_tip_state->GetCellCount(), 2u);

        // archiving test 
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "cell_mutation.arch";
        {
            AbstractCellProperty* const p_branching_state = new BranchingCellMutationState();
            p_branching_state->IncrementCellCount();
            AbstractCellProperty* const p_branching_state = new VesselCellMutationState();
            p_vessel_state->IncrementCellCount();
            AbstractCellProperty* const p_branching_state = new TipCellMutationState();
            p_tip_state->IncrementCellCount();

            TS_ASSERT_EQUALS(p_branching_state->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(p_vessel_state->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(p_tip_state->GetCellCount(), 1u);

            TS_ASSERT_EQUALS(dynamic_cast<AbstractCellMutationState*>(p_branching_state)->GetColour(), 2u);
            TS_ASSERT_EQUALS(dynamic_cast<AbstractCellMutationState*>(p_vessel_state)->GetColour(), 1u);
            TS_ASSERT_EQUALS(dynamic_cast<AbstractCellMutationState*>(p_tip_state)->GetColour(), 0u);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_branching_state;
            output_arch << p_vessel_state;
            output_arch << p_tip_state;

            delete p_branching_state;
            delete p_vessel_state;
            delete p_tip_state;
        }
        {
            AbstractCellProperty* p_arch_branching_state;
            AbstractCellProperty* p_arch_vessel_state;
            AbstractCellProperty* p_arch_tip_state;

            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_arch_branching_state;
            input_arch >> p_arch_vessel_state;
            input_arch >> p_arch_tip_state;

            TS_ASSERT_EQUALS(p_arch_branching_state->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(p_arch_vessel_state->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(p_arch_tip_state->GetCellCount(), 1u);

            BranchingCellMutationState* p_real_branching_state = dynamic_cast<BranchingCellMutationState*>(p_arch_branching_state);
            TS_ASSERT(p_real_branching_state != NULL);
            TS_ASSERT_EQUALS(p_real_branching_state->GetColour(), 2u);
            VesselCellMutationState* p_real_vessel_state = dynamic_cast<VesselCellMutationState*>(p_arch_vessel_state);
            TS_ASSERT(p_real_vessel_state != NULL);
            TS_ASSERT_EQUALS(p_real_vessel_state->GetColour(), 1u);
            TipCellMutationState* p_real_tip_state = dynamic_cast<TipCellMutationState*>(p_arch_tip_state);
            TS_ASSERT(p_real_tip_state != NULL);
            TS_ASSERT_EQUALS(p_real_tip_state->GetColour(), 0u);

            delete p_arch_branching_state;
            delete p_arch_vessel_state;
            delete p_arch_tip_state;
        }
    }

    void TestPinnedCellBoundaryCondition()
    {
        // creation of the mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0u, false, 10, 5));
        nodes.push_back(new Node<2>(1u, false, 9, 5));
        nodes.push_back(new Node<2>(2u, false, 8, 5));

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

        // fully constrain the first cell using the boundary condition 
        unsigned node_index_tip_cell = cell_population.GetLocationIndexUsingCell(0);
        std::vector<unsigned> pinned_node_indices;
        pinned_node_indices.push_back(node_index_tip_cell);
        typedef PinnedCellsBoundaryCondition<2,2> PinnedCellsBoundaryCondition;
        MAKE_PTR_ARGS(PinnedCellsBoundaryCondition, p_boundary_condition, (&cell_population, pinned_node_indices));

        // we check that the 
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