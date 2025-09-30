#ifndef TESTFORCESMODEL_HPP_
#define TESTFORCESMODEL_HPP_

#include <cxxtest/TestSuite.h>

#include "CellPropertyRegistry.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <iostream>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include "ProcessSpecificArchive.hpp"
#include "OutputFileHandler.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "ArchiveOpener.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "Debug.hpp"
#include "CheckpointArchiveTypes.hpp"

#include "AbstractCellBasedSimulationModifier.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "AbstractCellMutationState.hpp"
#include "AbstractCellPopulationBoundaryCondition.hpp"
#include "AbstractBoxDomainPdeModifier.hpp"
#include "AnastomosisWriter.hpp"
#include "AveragedSourceParabolicPde.hpp"
#include "Cell.hpp"
#include "CellBasedEventHandler.hpp"
#include "CellData.hpp"
#include "CellId.hpp"
#include "CellLabel.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FileComparison.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "ParabolicBoxDomainPdeModifier.hpp"
#include "ReplicatableVector.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "SmartPointers.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include "UniformSourceParabolicPde.hpp"

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

#include "BranchingSegmentMutationState.hpp"
#include "VesselTipMutationState.hpp"
#include "VesselSegmentMutationState.hpp"

#include "VegfEquationPde.hpp"
#include "VegfBoundaryCondition.hpp"
#include "MolecularConcentrationsDomainPdeModifier.hpp"

class TestForcesModel : public AbstractCellBasedTestSuite
{

public:

    void NoTestCellMutationStates()
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

        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        cells_generator.GenerateBasicRandom(cells, 5, p_differentiated_type);

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);
        MAKE_PTR(VesselTipMutationState, p_tip_state);

        // test the pointers: colours and type 
        TS_ASSERT_EQUALS(p_branching_state->GetColour(), 2u);
        TS_ASSERT_EQUALS(p_branching_state->IsType<VesselTipMutationState>(), false);
        TS_ASSERT_EQUALS(p_branching_state->IsType<VesselSegmentMutationState>(), false);
        TS_ASSERT_EQUALS(p_branching_state->IsType<BranchingSegmentMutationState>(), true);

        TS_ASSERT_EQUALS(p_vessel_state->GetColour(), 1u);
        TS_ASSERT_EQUALS(p_vessel_state->IsType<VesselTipMutationState>(), false);
        TS_ASSERT_EQUALS(p_vessel_state->IsType<VesselSegmentMutationState>(), true);
        TS_ASSERT_EQUALS(p_vessel_state->IsType<BranchingSegmentMutationState>(), false);

        TS_ASSERT_EQUALS(p_tip_state->GetColour(), 0u);
        TS_ASSERT_EQUALS(p_tip_state->IsType<VesselTipMutationState>(), true);
        TS_ASSERT_EQUALS(p_tip_state->IsType<VesselSegmentMutationState>(), false);
        TS_ASSERT_EQUALS(p_tip_state->IsType<BranchingSegmentMutationState>(), false);

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
            AbstractCellProperty* const p_branching_state = new BranchingSegmentMutationState();
            p_branching_state->IncrementCellCount();
            AbstractCellProperty* const p_vessel_state = new VesselSegmentMutationState();
            p_vessel_state->IncrementCellCount();
            AbstractCellProperty* const p_tip_state = new VesselTipMutationState();
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

            BranchingSegmentMutationState* p_real_branching_state = dynamic_cast<BranchingSegmentMutationState*>(p_arch_branching_state);
            TS_ASSERT(p_real_branching_state != NULL);
            TS_ASSERT_EQUALS(p_real_branching_state->GetColour(), 2u);
            VesselSegmentMutationState* p_real_vessel_state = dynamic_cast<VesselSegmentMutationState*>(p_arch_vessel_state);
            TS_ASSERT(p_real_vessel_state != NULL);
            TS_ASSERT_EQUALS(p_real_vessel_state->GetColour(), 1u);
            VesselTipMutationState* p_real_tip_state = dynamic_cast<VesselTipMutationState*>(p_arch_tip_state);
            TS_ASSERT(p_real_tip_state != NULL);
            TS_ASSERT_EQUALS(p_real_tip_state->GetColour(), 0u);

            delete p_arch_branching_state;
            delete p_arch_vessel_state;
            delete p_arch_tip_state;
        }
    }

    void NoTestVegfEquationPdeIn1D()
    {
        // creation of the mesh
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0u, false, 10));
        nodes.push_back(new Node<1>(1u, false, 9));
        nodes.push_back(new Node<1>(2u, false, 8));

        NodesOnlyMesh<1> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<1> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        // Chaste points corresponding to the vessel tip and a random point 
        ChastePoint<1> lower(0.0);
        ChastePoint<1> vessel_tip(8.0);
        double u = 1.0;

        VegfEquationPde<1> vegfpde(cell_population, 1.0, 1.0, 1.0, 0.1, 0.01);

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> fe_mesh;
        fe_mesh.ConstructFromMeshReader(mesh_reader);
        fe_mesh.Scale(10.0);
        std::map<CellPtr, unsigned> cell_pde_element_map;
        for (AbstractCellPopulation<1>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_pde_element_map[*cell_iter] = 0;
            if((*cell_iter)->GetMutationState()->IsType<VesselTipMutationState>()){
                cell_pde_element_map[*cell_iter] = 1;
            }
        }
        vegfpde.SetupSourceTerms(fe_mesh, &cell_pde_element_map);

        TS_ASSERT_DELTA(vegfpde.ComputeSourceTerm(lower, u, fe_mesh.GetElement(0)), -0.9, 1e-4);
        TS_ASSERT_DELTA(vegfpde.ComputeSourceTerm(vessel_tip, u, fe_mesh.GetElement(1)), -0.91, 1e-12);

        TS_ASSERT_DELTA(vegfpde.ComputeDuDtCoefficientFunction(lower), 1.0, 1e-12);
        TS_ASSERT_DELTA(vegfpde.ComputeDuDtCoefficientFunction(vessel_tip), 1.0, 1e-12);

        c_matrix<double, 1, 1> diffusion_lower = vegfpde.ComputeDiffusionTerm(lower, fe_mesh.GetElement(0));
        c_matrix<double, 1, 1> diffusion_vesseltip = vegfpde.ComputeDiffusionTerm(vessel_tip, fe_mesh.GetElement(1));

        TS_ASSERT_DELTA(diffusion_lower(0,0), 1.0, 1e-12);
        TS_ASSERT_DELTA(diffusion_vesseltip(0,0), 1.0, 1e-12);

        VegfBoundaryCondition<1> vegfboundarycondition(0.1, 8);
        TS_ASSERT_DELTA(vegfboundarycondition.GetValue(lower), 0.0, 1e-3);
        TS_ASSERT_DELTA(vegfboundarycondition.GetValue(vessel_tip), 0.1, 1e-3);

        // archive test
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "vegfequationpde.arch";
        {
            // Create a PDE object and the boundary condition 
            AbstractLinearParabolicPde<1,1>* const p_vegfpde = new VegfEquationPde<1>(cell_population, 1.0, 1.0, 1.0, 0.1, 0.01);

            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* output_arch = arch_opener.GetCommonArchive();
            (*output_arch) << p_vegfpde;

            delete p_vegfpde;
        }
        {
            AbstractLinearParabolicPde<1,1>* p_vegfpde;

            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* input_arch = arch_opener.GetCommonArchive();
            (*input_arch) >> p_vegfpde;

            // Test that the PDE and its member variables were archived correctly
            TS_ASSERT(static_cast<VegfEquationPde<1>*>(p_vegfpde) != NULL);

            VegfEquationPde<1>* p_static_cast_vegfpde = static_cast<VegfEquationPde<1>*>(p_vegfpde);
            TS_ASSERT_EQUALS(p_static_cast_vegfpde->GetDuDtCoefficient(),1.0);
            TS_ASSERT_EQUALS(p_static_cast_vegfpde->GetDiffusionCoefficient(),1.0);
            TS_ASSERT_EQUALS(p_static_cast_vegfpde->GetDecayCoefficient(),1.0);
            TS_ASSERT_EQUALS(p_static_cast_vegfpde->GetCreationCoefficient(),0.1);
            TS_ASSERT_EQUALS(p_static_cast_vegfpde->GetConsumptionCoefficient(),0.01);

            // Avoid memory leaks
            delete p_vegfpde;
        }
    }

    void NoTestVegfEquationPdeIn2D()
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
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

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

        // Chaste points corresponding to the vessel tip and a random point 
        ChastePoint<2> lower(0.0, 0.0);
        ChastePoint<2> vessel_tip(8.0, 1.0);
        double u = 1.0;

        VegfEquationPde<2> vegfpde(cell_population, 1.0, 1.0, 1.0, 0.1, 0.01);

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        TetrahedralMesh<2,2> fe_mesh;
        fe_mesh.ConstructFromMeshReader(mesh_reader);
        fe_mesh.Scale(2.0, 1.0);
        std::map<CellPtr, unsigned> cell_pde_element_map;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_pde_element_map[*cell_iter] = 0;
            if((*cell_iter)->GetMutationState()->IsType<VesselTipMutationState>()){
                cell_pde_element_map[*cell_iter] = 1;
            }
        }
        vegfpde.SetupSourceTerms(fe_mesh, &cell_pde_element_map);

        TS_ASSERT_DELTA(vegfpde.ComputeSourceTerm(lower, u, fe_mesh.GetElement(0)), -0.9, 1e-3);
        TS_ASSERT_DELTA(vegfpde.ComputeSourceTerm(vessel_tip, u, fe_mesh.GetElement(1)), -0.91, 1e-3);

        TS_ASSERT_DELTA(vegfpde.ComputeDuDtCoefficientFunction(lower), 1.0, 1e-12);
        TS_ASSERT_DELTA(vegfpde.ComputeDuDtCoefficientFunction(vessel_tip), 1.0, 1e-12);

        c_matrix<double, 2, 2> diffusion_lower = vegfpde.ComputeDiffusionTerm(lower, fe_mesh.GetElement(0));
        c_matrix<double, 2, 2> diffusion_vesseltip = vegfpde.ComputeDiffusionTerm(vessel_tip, fe_mesh.GetElement(1));

        TS_ASSERT_DELTA(diffusion_lower(0,0), 1.0, 1e-12);
        TS_ASSERT_DELTA(diffusion_lower(1,1), 1.0, 1e-12);
        TS_ASSERT_DELTA(diffusion_lower(0,1), 0, 1e-12);
        TS_ASSERT_DELTA(diffusion_lower(1,0), 0, 1e-12);
        TS_ASSERT_DELTA(diffusion_vesseltip(0,0), 1.0, 1e-12);
        TS_ASSERT_DELTA(diffusion_vesseltip(1,1), 1.0, 1e-12);
        TS_ASSERT_DELTA(diffusion_vesseltip(0,1), 0, 1e-12);
        TS_ASSERT_DELTA(diffusion_vesseltip(1,0), 0, 1e-12);

        VegfBoundaryCondition<2> vegfboundarycondition(0.1, 8);
        TS_ASSERT_DELTA(vegfboundarycondition.GetValue(lower), 0.0, 1e-3);
        TS_ASSERT_DELTA(vegfboundarycondition.GetValue(vessel_tip), 0.1, 1e-3);

        // archive test
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "vegfequationpde.arch";
        {
            // Create a PDE object and the boundary condition 
            AbstractLinearParabolicPde<2,2>* const p_vegfpde = new VegfEquationPde<2>(cell_population, 1.0, 1.0, 1.0, 0.1, 0.01);

            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* output_arch = arch_opener.GetCommonArchive();
            (*output_arch) << p_vegfpde;

            delete p_vegfpde;
        }
        {
            AbstractLinearParabolicPde<2,2>* p_vegfpde;

            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* input_arch = arch_opener.GetCommonArchive();
            (*input_arch) >> p_vegfpde;

            // Test that the PDE and its member variables were archived correctly
            TS_ASSERT(static_cast<VegfEquationPde<2>*>(p_vegfpde) != NULL);

            VegfEquationPde<2>* p_static_cast_vegfpde = static_cast<VegfEquationPde<2>*>(p_vegfpde);
            TS_ASSERT_EQUALS(p_static_cast_vegfpde->GetDuDtCoefficient(),1.0);
            TS_ASSERT_EQUALS(p_static_cast_vegfpde->GetDiffusionCoefficient(),1.0);
            TS_ASSERT_EQUALS(p_static_cast_vegfpde->GetDecayCoefficient(),1.0);
            TS_ASSERT_EQUALS(p_static_cast_vegfpde->GetCreationCoefficient(),0.1);
            TS_ASSERT_EQUALS(p_static_cast_vegfpde->GetConsumptionCoefficient(),0.01);

            // Avoid memory leaks
            delete p_vegfpde;
        }
    }

    void NoTestVegfEquationPdeIn3D()
    {
        // creation of the mesh
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u, false, 10, 5, 5));
        nodes.push_back(new Node<3>(1u, false, 9, 5, 5));
        nodes.push_back(new Node<3>(2u, false, 8, 5, 5));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        // Chaste points corresponding to the vessel tip and a random point 
        ChastePoint<3> lower(0.0, 0.0, 0.0);
        ChastePoint<3> vessel_tip(8.0, 0.0, 0.0);
        double u = 1.0;

        VegfEquationPde<3> vegfpde(cell_population, 1.0, 1.0, 1.0, 0.1, 0.01);

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        TetrahedralMesh<3,3> fe_mesh;
        fe_mesh.ConstructFromMeshReader(mesh_reader);
        fe_mesh.Scale(136.0, 1.0, 1.0);
        std::map<CellPtr, unsigned> cell_pde_element_map;
        for (AbstractCellPopulation<3>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_pde_element_map[*cell_iter] = 0;
            if((*cell_iter)->GetMutationState()->IsType<VesselTipMutationState>()){
                cell_pde_element_map[*cell_iter] = 1;
            }
        }
        vegfpde.SetupSourceTerms(fe_mesh, &cell_pde_element_map);

        TS_ASSERT_DELTA(vegfpde.ComputeSourceTerm(lower, u, fe_mesh.GetElement(0)), -0.9, 1e-4);
        TS_ASSERT_DELTA(vegfpde.ComputeSourceTerm(vessel_tip, u, fe_mesh.GetElement(8)), -0.91, 1e-2);

        TS_ASSERT_DELTA(vegfpde.ComputeDuDtCoefficientFunction(lower), 1.0, 1e-12);
        TS_ASSERT_DELTA(vegfpde.ComputeDuDtCoefficientFunction(vessel_tip), 1.0, 1e-12);

        c_matrix<double, 3, 3> diffusion_lower = vegfpde.ComputeDiffusionTerm(lower, fe_mesh.GetElement(0));
        c_matrix<double, 3, 3> diffusion_vesseltip = vegfpde.ComputeDiffusionTerm(vessel_tip, fe_mesh.GetElement(1));

        TS_ASSERT_DELTA(diffusion_lower(0,0), 1.0, 1e-12);
        TS_ASSERT_DELTA(diffusion_lower(1,1), 1.0, 1e-12);
        TS_ASSERT_DELTA(diffusion_lower(2,2), 1.0, 1e-12);
        TS_ASSERT_DELTA(diffusion_lower(0,1), 0, 1e-12);
        TS_ASSERT_DELTA(diffusion_lower(0,2), 0, 1e-12);
        TS_ASSERT_DELTA(diffusion_lower(1,0), 0, 1e-12);
        TS_ASSERT_DELTA(diffusion_lower(1,2), 0, 1e-12);
        TS_ASSERT_DELTA(diffusion_lower(2,0), 0, 1e-12);
        TS_ASSERT_DELTA(diffusion_lower(2,1), 0, 1e-12);
        TS_ASSERT_DELTA(diffusion_vesseltip(0,0), 1.0, 1e-12);
        TS_ASSERT_DELTA(diffusion_vesseltip(1,1), 1.0, 1e-12);
        TS_ASSERT_DELTA(diffusion_vesseltip(2,2), 1.0, 1e-12);
        TS_ASSERT_DELTA(diffusion_vesseltip(0,1), 0, 1e-12);
        TS_ASSERT_DELTA(diffusion_vesseltip(0,2), 0, 1e-12);
        TS_ASSERT_DELTA(diffusion_vesseltip(1,0), 0, 1e-12);
        TS_ASSERT_DELTA(diffusion_vesseltip(1,2), 0, 1e-12);
        TS_ASSERT_DELTA(diffusion_vesseltip(2,0), 0, 1e-12);
        TS_ASSERT_DELTA(diffusion_vesseltip(2,1), 0, 1e-12);

        VegfBoundaryCondition<3> vegfboundarycondition(0.1, 8);
        TS_ASSERT_DELTA(vegfboundarycondition.GetValue(lower), 0.0, 1e-3);
        TS_ASSERT_DELTA(vegfboundarycondition.GetValue(vessel_tip), 0.1, 1e-3);

        // archive test
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "vegfequationpde.arch";
        {
            // Create a PDE object and the boundary condition 
            AbstractLinearParabolicPde<3,3>* const p_vegfpde = new VegfEquationPde<3>(cell_population, 1.0, 1.0, 1.0, 0.1, 0.01);

            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* output_arch = arch_opener.GetCommonArchive();
            (*output_arch) << p_vegfpde;

            delete p_vegfpde;
        }
        {
            AbstractLinearParabolicPde<3,3>* p_vegfpde;

            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* input_arch = arch_opener.GetCommonArchive();
            (*input_arch) >> p_vegfpde;

            // Test that the PDE and its member variables were archived correctly
            TS_ASSERT(static_cast<VegfEquationPde<3>*>(p_vegfpde) != NULL);

            VegfEquationPde<3>* p_static_cast_vegfpde = static_cast<VegfEquationPde<3>*>(p_vegfpde);
            TS_ASSERT_EQUALS(p_static_cast_vegfpde->GetDuDtCoefficient(),1.0);
            TS_ASSERT_EQUALS(p_static_cast_vegfpde->GetDiffusionCoefficient(),1.0);
            TS_ASSERT_EQUALS(p_static_cast_vegfpde->GetDecayCoefficient(),1.0);
            TS_ASSERT_EQUALS(p_static_cast_vegfpde->GetCreationCoefficient(),0.1);
            TS_ASSERT_EQUALS(p_static_cast_vegfpde->GetConsumptionCoefficient(),0.01);

            // Avoid memory leaks
            delete p_vegfpde;
        }
    }

    void NoTestMolecularConcentrationsDomainPdeModifierIn1D() throw(Exception)
    {
        // creation of the mesh
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0u, false, 10));
        nodes.push_back(new Node<1>(1u, false, 9));
        nodes.push_back(new Node<1>(2u, false, 8));

        NodesOnlyMesh<1> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<1> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        // Chaste points corresponding to the vessel tip and a random point 
        ChastePoint<1> lower(0.0);
        ChastePoint<1> upper(10.0);
        MAKE_PTR_ARGS(ChasteCuboid<1>, p_cuboid, (lower, upper));

        TS_ASSERT_DELTA(p_cuboid->rGetUpperCorner()[0],10,1e-5);
        TS_ASSERT_DELTA(p_cuboid->rGetLowerCorner()[0],0,1e-5);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(VegfEquationPde<1>, p_pde, (cell_population, 1.0, 1.0, 1.0, 0.1, 0.01));
        MAKE_PTR_ARGS(VegfBoundaryCondition<1>, p_bc, (0.1, 0.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(MolecularConcentrationsDomainPdeModifier<1>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 1.0, nullptr, 0.0, 1.0, 0.1));

        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumNodes(),11u);
        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumBoundaryNodes(),2u);
        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumElements(),10u);
        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumBoundaryElements(),2u);

        // Test that member variables are initialised correctly
        p_pde_modifier->SetDependentVariableName("vegf_femesh_variable");
        TS_ASSERT_EQUALS(p_pde_modifier->rGetDependentVariableName(), "vegf_femesh_variable");

        // Coverage
        TS_ASSERT_EQUALS(p_pde_modifier->GetOutputGradient(),false); // Defaults to false
        p_pde_modifier->SetOutputGradient(true);
        TS_ASSERT_EQUALS(p_pde_modifier->GetOutputGradient(),true);

        // setup initial vector
        p_pde_modifier->SetupInitialSolutionVector(cell_population);
        Vec solution = p_pde_modifier->GetSolution();
        ReplicatableVector solution_repl(solution);
        TS_ASSERT_EQUALS(solution_repl.GetSize(), 11u);
        for (unsigned i=0; i<11; i++)
        {
            if(p_pde_modifier->mpFeMesh->GetNode(i)->rGetLocation()[0] == 0.0){
                TS_ASSERT_DELTA(solution_repl[i], 1.0, 1e-6); // source term
            } else {
                TS_ASSERT_DELTA(solution_repl[i], 0.1, 1e-6); // background
            }
        }

        std::string archive_file = "molecularconcentrationsdomainpdemodifier.arch";
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + archive_file;
        {
            // Create PDE and boundary condition objects
            MAKE_PTR_ARGS(VegfEquationPde<1>, p_pde, (cell_population, 1.0, 1.0, 1.0, 0.1, 0.01));
            MAKE_PTR_ARGS(VegfBoundaryCondition<1>, p_bc, (0.1, 8));

            // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
            ChastePoint<1> lower(0.0);
            ChastePoint<1> upper(10.0);
            MAKE_PTR_ARGS(ChasteCuboid<1>, p_cuboid, (lower, upper));

            // Create a PDE modifier and set the name of the dependent variable in the PDE
            MolecularConcentrationsDomainPdeModifier<1> modifier(p_pde, p_bc, false, p_cuboid, 1.0, nullptr, 0.0, 1.0, 0.1);
            modifier.SetDependentVariableName("vegf_femesh_variable");
            modifier.SetupInitialSolutionVector(cell_population);
            modifier.SetOutputGradient(true);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            AbstractBoxDomainPdeModifier<1>* const p_modifier = &modifier;
            output_arch << p_modifier;
        }
        {
            AbstractBoxDomainPdeModifier<1>* p_modifier;

            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_modifier;

            TS_ASSERT_EQUALS((static_cast<MolecularConcentrationsDomainPdeModifier<1>*>(p_modifier))->rGetDependentVariableName(), "vegf_femesh_variable");
            TS_ASSERT_DELTA((static_cast<MolecularConcentrationsDomainPdeModifier<1>*>(p_modifier))->GetStepSize(), 1.0, 1e-5);
            TS_ASSERT_EQUALS((static_cast<MolecularConcentrationsDomainPdeModifier<1>*>(p_modifier))->AreBcsSetOnBoxBoundary(), true);

            // comparison solution vector 
            Vec solution = (static_cast<MolecularConcentrationsDomainPdeModifier<1>*>(p_modifier))->GetSolution();
            ReplicatableVector solution_repl(solution);
            TS_ASSERT_EQUALS(solution_repl.GetSize(), 11u);
            for (unsigned i=0; i<11; i++)
            {
                if((static_cast<MolecularConcentrationsDomainPdeModifier<1>*>(p_modifier))->mpFeMesh->GetNode(i)->rGetLocation()[0] == 0.0){
                    TS_ASSERT_DELTA(solution_repl[i], 1.0, 1e-6); // source term
                } else {
                    TS_ASSERT_DELTA(solution_repl[i], 0.1, 1e-6); // background
                }
            }

            delete p_modifier;
        }
    }

    void NoTestMolecularConcentrationsDomainPdeModifierIn2D() throw(Exception)
    {
        // creation of the mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0u, false, 10, 0));
        nodes.push_back(new Node<2>(1u, false, 9, 0));
        nodes.push_back(new Node<2>(2u, false, 8, 0));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

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

        // Chaste points corresponding to the vessel tip and a random point 
        ChastePoint<2> lower(0.0, 0.0);
        ChastePoint<2> upper(10.0, 10.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        TS_ASSERT_DELTA(p_cuboid->rGetUpperCorner()[0],10,1e-5);
        TS_ASSERT_DELTA(p_cuboid->rGetUpperCorner()[1],10,1e-5);
        TS_ASSERT_DELTA(p_cuboid->rGetLowerCorner()[0],0,1e-5);
        TS_ASSERT_DELTA(p_cuboid->rGetLowerCorner()[1],0,1e-5);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(VegfEquationPde<2>, p_pde, (cell_population, 1.0, 1.0, 1.0, 0.1, 0.01));
        MAKE_PTR_ARGS(VegfBoundaryCondition<2>, p_bc, (0.1, 0.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(MolecularConcentrationsDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 1.0,nullptr, 0.0, 1.0, 0.1));

        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumNodes(),121u);
        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumBoundaryNodes(),40u);
        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumElements(),200u);
        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumBoundaryElements(),40u);

        // Test that member variables are initialised correctly
        p_pde_modifier->SetDependentVariableName("vegf_femesh_variable");
        TS_ASSERT_EQUALS(p_pde_modifier->rGetDependentVariableName(), "vegf_femesh_variable");

        // Coverage
        TS_ASSERT_EQUALS(p_pde_modifier->GetOutputGradient(),false); // Defaults to false
        p_pde_modifier->SetOutputGradient(true);
        TS_ASSERT_EQUALS(p_pde_modifier->GetOutputGradient(),true);

        // setup initial vector
        p_pde_modifier->SetupInitialSolutionVector(cell_population);
        Vec solution = p_pde_modifier->GetSolution();
        ReplicatableVector solution_repl(solution);
        TS_ASSERT_EQUALS(solution_repl.GetSize(), 121u);
        for (unsigned i=0; i<121; i++)
        {
            if(p_pde_modifier->mpFeMesh->GetNode(i)->rGetLocation()[0] == 0.0){
                TS_ASSERT_DELTA(solution_repl[i], 1.0, 1e-6); // source term
            } else {
                TS_ASSERT_DELTA(solution_repl[i], 0.1, 1e-6); // background
            }
        }

        std::string archive_file = "molecularconcentrationsdomainpdemodifier.arch";
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + archive_file;
        // {
        //     // Create PDE and boundary condition objects
        //     MAKE_PTR_ARGS(VegfEquationPde<2>, p_pde, (cell_population, 1.0, 1.0, 1.0, 0.1, 0.01));
        //     MAKE_PTR_ARGS(VegfBoundaryCondition<2>, p_bc, (0.1, 8));

        //     // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        //     ChastePoint<2> lower(0.0, 0.0);
        //     ChastePoint<2> upper(10.0, 10.0);
        //     MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        //     // Create a PDE modifier and set the name of the dependent variable in the PDE
        //     MolecularConcentrationsDomainPdeModifier<2> modifier(p_pde, p_bc, false, p_cuboid, 1.0, nullptr, 0.0,0.1, 0.1);
        //     modifier.SetDependentVariableName("vegf_femesh_variable");

        //     ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
        //     boost::archive::text_oarchive* output_arch = arch_opener.GetCommonArchive();

        //     // Serialize via pointer
        //     AbstractCellBasedSimulationModifier<2,2>* const p_modifier = &modifier;
        //     (*output_arch) << p_modifier;
        // }
        // {
        //     AbstractCellBasedSimulationModifier<2,2>* p_modifier2;

        //     ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
        //     boost::archive::text_iarchive* input_arch = arch_opener.GetCommonArchive();

        //     (*input_arch) >> p_modifier2;
        //     TRACE("Test if we can access archive")

        //     delete p_modifier2;
        // }
    }

    void NoTestMolecularConcentrationsDomainPdeModifierIn3D() throw(Exception)
    {
        // creation of the mesh
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u, false, 10, 0));
        nodes.push_back(new Node<3>(1u, false, 9, 0));
        nodes.push_back(new Node<3>(2u, false, 8, 0));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        // Chaste points corresponding to the vessel tip and a random point 
        ChastePoint<3> lower(0.0, 0.0, 0.0);
        ChastePoint<3> upper(10.0, 10.0, 10.0);
        MAKE_PTR_ARGS(ChasteCuboid<3>, p_cuboid, (lower, upper));

        TS_ASSERT_DELTA(p_cuboid->rGetUpperCorner()[0],10,1e-5);
        TS_ASSERT_DELTA(p_cuboid->rGetUpperCorner()[1],10,1e-5);
        TS_ASSERT_DELTA(p_cuboid->rGetUpperCorner()[2],10,1e-5);
        TS_ASSERT_DELTA(p_cuboid->rGetLowerCorner()[0],0,1e-5);
        TS_ASSERT_DELTA(p_cuboid->rGetLowerCorner()[1],0,1e-5);
        TS_ASSERT_DELTA(p_cuboid->rGetLowerCorner()[2],0,1e-5);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(VegfEquationPde<3>, p_pde, (cell_population, 1.0, 1.0, 1.0, 0.1, 0.01));
        MAKE_PTR_ARGS(VegfBoundaryCondition<3>, p_bc, (0.1, 0.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(MolecularConcentrationsDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 1.0,nullptr, 0.0, 1.0, 0.1));

        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumNodes(),1331u);
        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumBoundaryNodes(),602u);
        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumElements(),6000u);
        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumBoundaryElements(),1200u);

        // Test that member variables are initialised correctly
        p_pde_modifier->SetDependentVariableName("vegf_femesh_variable");
        TS_ASSERT_EQUALS(p_pde_modifier->rGetDependentVariableName(), "vegf_femesh_variable");

        // Coverage
        TS_ASSERT_EQUALS(p_pde_modifier->GetOutputGradient(),false); // Defaults to false
        p_pde_modifier->SetOutputGradient(true);
        TS_ASSERT_EQUALS(p_pde_modifier->GetOutputGradient(),true);

        // setup initial vector
        p_pde_modifier->SetupInitialSolutionVector(cell_population);
        Vec solution = p_pde_modifier->GetSolution();
        ReplicatableVector solution_repl(solution);
        TS_ASSERT_EQUALS(solution_repl.GetSize(), 1331u);
        for (unsigned i=0; i<1331; i++)
        {
            if(p_pde_modifier->mpFeMesh->GetNode(i)->rGetLocation()[0] == 0.0){
                TS_ASSERT_DELTA(solution_repl[i], 1.0, 1e-6); // source term
            } else {
                TS_ASSERT_DELTA(solution_repl[i], 0.1, 1e-6); // background
            }
        }

        std::string archive_file = "molecularconcentrationsdomainpdemodifier.arch";
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + archive_file;
        // {
        //     // Create PDE and boundary condition objects
        //     MAKE_PTR_ARGS(VegfEquationPde<3>, p_pde, (cell_population, 1.0, 1.0, 1.0, 0.1, 0.01));
        //     MAKE_PTR_ARGS(VegfBoundaryCondition<3>, p_bc, (0.1, 8));

        //     // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        //     ChastePoint<3> lower(0.0, 0.0, 0.0);
        //     ChastePoint<3> upper(10.0, 10.0, 10.0);
        //     MAKE_PTR_ARGS(ChasteCuboid<3>, p_cuboid, (lower, upper));

        //     // Create a PDE modifier and set the name of the dependent variable in the PDE
        //     MolecularConcentrationsDomainPdeModifier<3> modifier(p_pde, p_bc, false, p_cuboid, 1.0, nullptr, 0.0,0.1, 0.1);
        //     modifier.SetDependentVariableName("vegf_femesh_variable");

        //     ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
        //     boost::archive::text_oarchive* output_arch = arch_opener.GetCommonArchive();

        //     // Serialize via pointer
        //     AbstractCellBasedSimulationModifier<3,3>* const p_modifier = &modifier;
        //     (*output_arch) << p_modifier;
        // }
        // {
        //     AbstractCellBasedSimulationModifier<3,3>* p_modifier2;

        //     ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
        //     boost::archive::text_iarchive* input_arch = arch_opener.GetCommonArchive();

        //     (*input_arch) >> p_modifier2;
        //     TRACE("Test if we can access archive")

        //     delete p_modifier2;
        // }
    }

    void NoTestRandomForceIn1D()
    {
        // creation of the mesh
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0u, false, 10));
        nodes.push_back(new Node<1>(1u, false, 9));
        nodes.push_back(new Node<1>(2u, false, 8));

        NodesOnlyMesh<1> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<1> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        // we create the random force object 
        RandomForce<1> random_force(0.1);

        random_force.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            TS_ASSERT_DELTA(abs(cell_population.GetNode(node_index)->rGetAppliedForce()[0]), 0.1, 1e-1);
            //TS_ASSERT_EQUALS(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.1); // random norm
            TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.1, 1e-5);
        }

        // test archiving class 
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "random_force.arch";
        {
            AbstractForce<1>* const p_force = new RandomForce<1>(0.1);
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_force;
            delete p_force;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractForce<1>* p_force;
            input_arch >> p_force;

            TS_ASSERT_EQUALS(dynamic_cast<RandomForce<1>*>(p_force)->GetRandomSensitivity(), 0.1);

            delete p_force;
        }
    }

    void NoTestRandomForceIn2D()
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
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

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
            TS_ASSERT_DELTA(abs(cell_population.GetNode(node_index)->rGetAppliedForce()[0]), 0.1, 1e-1);
            TS_ASSERT_DELTA(abs(cell_population.GetNode(node_index)->rGetAppliedForce()[1]), 0.1, 1e-1);
            //TS_ASSERT_EQUALS(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.1); // random norm
            TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.1, 1e-5);
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

    void NoTestRandomForceIn3D()
    {
        // creation of the mesh
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u, false, 10, 5, 5));
        nodes.push_back(new Node<3>(1u, false, 9, 5, 5));
        nodes.push_back(new Node<3>(2u, false, 8, 5, 5));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        // we create the random force object 
        RandomForce<3> random_force(0.1);

        random_force.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            TS_ASSERT_DELTA(abs(cell_population.GetNode(node_index)->rGetAppliedForce()[0]), 0.1, 1e-1);
            TS_ASSERT_DELTA(abs(cell_population.GetNode(node_index)->rGetAppliedForce()[1]), 0.1, 1e-1);
            TS_ASSERT_DELTA(abs(cell_population.GetNode(node_index)->rGetAppliedForce()[2]), 0.1, 1e-1);
            //TS_ASSERT_EQUALS(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.1); // random norm
            TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.1, 1e-5);
        }

        // test archiving class 
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "random_force.arch";
        {
            AbstractForce<3>* const p_force = new RandomForce<3>(0.1);
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_force;
            delete p_force;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractForce<3>* p_force;
            input_arch >> p_force;

            TS_ASSERT_EQUALS(dynamic_cast<RandomForce<3>*>(p_force)->GetRandomSensitivity(), 0.1);

            delete p_force;
        }
    }

    void NoTestMechanicalForceIn1D()
    {
        // creation of the mesh
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0u, false, 10));
        nodes.push_back(new Node<1>(1u, false, 9));
        nodes.push_back(new Node<1>(2u, false, 8));

        NodesOnlyMesh<1> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<1> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        // we create the random force object 
        LinearMechanicalForceModified<1> mechanical_force;

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
            AbstractForce<1>* const p_force = new LinearMechanicalForceModified<1>();
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_force;
            delete p_force;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractForce<1>* p_force;
            input_arch >> p_force;

            delete p_force;
        }
    }

    void NoTestMechanicalForceIn2D()
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
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

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

    void NoTestMechanicalForceIn3D()
    {
        // creation of the mesh
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u, false, 10, 5, 5));
        nodes.push_back(new Node<3>(1u, false, 9, 5, 5));
        nodes.push_back(new Node<3>(2u, false, 8, 5, 5));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        // we create the random force object 
        LinearMechanicalForceModified<3> mechanical_force;

        mechanical_force.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 1e-1, 1e-1); 
        }

        // test archiving class 
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "random_force.arch";
        {
            AbstractForce<3>* const p_force = new LinearMechanicalForceModified<3>();
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_force;
            delete p_force;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractForce<3>* p_force;
            input_arch >> p_force;

            delete p_force;
        }
    }

    void NoTestChemotacticForceIn1D()
    {
        // creation of the mesh
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0u, false, 10));
        nodes.push_back(new Node<1>(1u, false, 9));
        nodes.push_back(new Node<1>(2u, false, 8));

        NodesOnlyMesh<1> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<1> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        // we create the random force object 
        ChemoForce<1> chemo_force(1e-4, 1e-4, 1e-4, 1e-4, 0.1, -0.1, 0.0);

        chemo_force.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            if(p_cell->GetMutationState()->IsType<VesselTipMutationState>()){
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[0], -1e-5, 1e-5);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 1e-5, 1e-5);
            } else {
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 0.0);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.0, 1.0);
            }
        }

        // test archiving class 
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "chemo_force.arch";
        {
            AbstractForce<1>* const p_force = new ChemoForce<1>(1e-4);
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_force;
            delete p_force;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractForce<1>* p_force;
            input_arch >> p_force;

            TS_ASSERT_EQUALS(dynamic_cast<ChemoForce<1>*>(p_force)->GetChemotacticSensitivity(), 0.1);

            delete p_force;
        }
    }

    void NoTestChemotacticForceIn2D()
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
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

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
        ChemoForce<2> chemo_force(1e-4, 1e-4, 1e-4, 1e-4, 0.1, -0.1, 0.0);

        chemo_force.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            if(p_cell->GetMutationState()->IsType<VesselTipMutationState>()){
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[0], -1e-5, 1e-5);
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 1e-5, 1e-5);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 1e-5, 1e-5); // random norm
            } else {
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 0.0);
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 0.0);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.0, 1.0); // random norm
            }
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

            TS_ASSERT_EQUALS(dynamic_cast<ChemoForce<2>*>(p_force)->GetChemotacticSensitivity(), 0.1);

            delete p_force;
        }
    }

    void NoTestChemotacticForceIn3D()
    {
        // creation of the mesh
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u, false, 10, 5, 5));
        nodes.push_back(new Node<3>(1u, false, 9, 5, 5));
        nodes.push_back(new Node<3>(2u, false, 8, 5, 5));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        // we create the random force object 
        ChemoForce<3> chemo_force(1e-4, 1e-4, 1e-4, 1e-4, 0.1, -0.1, -0.1);

        chemo_force.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            if(p_cell->GetMutationState()->IsType<VesselTipMutationState>()){
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[0], -1e-5, 1e-5);
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 1e-5, 1e-5);
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[2], 1e-5, 1e-5);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 1e-5, 1e-5);
            } else {
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 0.0);
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 0.0);
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[2], 0.0);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.0, 1.0);
            }
        }

        // test archiving class 
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "chemo_force.arch";
        {
            AbstractForce<3>* const p_force = new ChemoForce<3>(1e-4);
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_force;
            delete p_force;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractForce<3>* p_force;
            input_arch >> p_force;

            TS_ASSERT_EQUALS(dynamic_cast<ChemoForce<3>*>(p_force)->GetChemotacticSensitivity(), 0.1);

            delete p_force;
        }
    }

    void NoTestChemotacticForceWithConstantVegfIn1D()
    {
        // creation of the mesh
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0u, false, 10));
        nodes.push_back(new Node<1>(1u, false, 9));
        nodes.push_back(new Node<1>(2u, false, 8));

        NodesOnlyMesh<1> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<1> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        // we create the random force object 
        ChemoForceWithConstantVegf<1> chemo_force_withconstantvegf(1e-4, 1e-4, 0.1);

        chemo_force_withconstantvegf.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            if(p_cell->GetMutationState()->IsType<VesselTipMutationState>()){
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[0], -1e-5, 1e-5);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 1e-5, 1e-5);
            } else {
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 0.0);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.0, 1.0); 
            }
        }

        // test archiving class 
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "chemo_force_withconstantvegf.arch";
        {
            AbstractForce<1>* const p_force = new ChemoForceWithConstantVegf<1>(1e-4);
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_force;
            delete p_force;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractForce<1>* p_force;
            input_arch >> p_force;

            TS_ASSERT_EQUALS(dynamic_cast<ChemoForceWithConstantVegf<1>*>(p_force)->GetChemotacticSensitivity(), 0.1);

            delete p_force;
        }
    }

    void NoTestChemotacticForceWithConstantVegfIn2D()
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
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

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
        ChemoForceWithConstantVegf<2> chemo_force_withconstantvegf(1e-4, 1e-4, 0.1);

        chemo_force_withconstantvegf.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            if(p_cell->GetMutationState()->IsType<VesselTipMutationState>()){
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[0], -1e-5, 1e-5);
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 0.0);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 1e-5, 1e-5); // random norm
            } else {
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 0.0);
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 0.0);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.0, 1.0); // random norm
            }
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

            TS_ASSERT_EQUALS(dynamic_cast<ChemoForceWithConstantVegf<2>*>(p_force)->GetChemotacticSensitivity(), 0.1);

            delete p_force;
        }
    }

    void NoTestChemotacticForceWithConstantVegfIn3D()
    {
        // creation of the mesh
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u, false, 10, 5, 5));
        nodes.push_back(new Node<3>(1u, false, 9, 5, 5));
        nodes.push_back(new Node<3>(2u, false, 8, 5, 5));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        // we create the random force object 
        ChemoForceWithConstantVegf<3> chemo_force_withconstantvegf(1e-4, 1e-4, 0.1);

        chemo_force_withconstantvegf.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            if(p_cell->GetMutationState()->IsType<VesselTipMutationState>()){
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[0], -1e-5, 1e-5);
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 0.0);
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[2], 0.0);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 1e-5, 1e-5);
            } else {
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 0.0);
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 0.0);
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[2], 0.0);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.0, 1.0); 
            }
        }

        // test archiving class 
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "chemo_force_withconstantvegf.arch";
        {
            AbstractForce<3>* const p_force = new ChemoForceWithConstantVegf<3>(1e-4);
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_force;
            delete p_force;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractForce<3>* p_force;
            input_arch >> p_force;

            TS_ASSERT_EQUALS(dynamic_cast<ChemoForceWithConstantVegf<3>*>(p_force)->GetChemotacticSensitivity(), 0.1);

            delete p_force;
        }
    }

    void NoTestChemotacticForceWithPdeIn1D() 
    {
        // creation of the mesh
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0u, false, 2, 5));
        nodes.push_back(new Node<1>(1u, false, 1, 5));
        nodes.push_back(new Node<1>(2u, false, 0, 5));

        NodesOnlyMesh<1> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<1> cell_population(mesh, cells);

        ChastePoint<1> lower(0.0);
        ChastePoint<1> upper(10.0);
        MAKE_PTR_ARGS(ChasteCuboid<1>, p_cuboid, (lower, upper));

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        MAKE_PTR_ARGS(VegfEquationPde<1>, p_pde, (cell_population, 1.0, 1.0, 1.0, 0.1, 0.01));
        MAKE_PTR_ARGS(VegfBoundaryCondition<1>, p_bc, (1.0, 0));
        MAKE_PTR_ARGS(MolecularConcentrationsDomainPdeModifier<1>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 1.0,nullptr,0.0,1.0, 0.1));
        p_pde_modifier->SetDependentVariableName("vegf_femesh_variable");
        p_pde_modifier->SetOutputGradient(true);
        p_pde_modifier->SetupInitialSolutionVector(cell_population);

        // we create the random force object 
        ChemoForceWithPdes<1> chemo_force_withpdes(0.1, 1e-4, 1e-4, 1e-4, p_pde_modifier);

        chemo_force_withpdes.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            if(p_cell->GetMutationState()->IsType<VesselTipMutationState>()){
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 1e-2, 1e-3);
            } else {
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.0, 1e-5); 
            }
        }

        // test archiving class 
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "chemo_force_withpdes.arch";
        {
            MAKE_PTR_ARGS(VegfEquationPde<1>, p_pde, (cell_population, 1.0, 1.0, 1.0, 0.1, 0.01));
            MAKE_PTR_ARGS(VegfBoundaryCondition<1>, p_bc, (1.0, 8));
            MAKE_PTR_ARGS(MolecularConcentrationsDomainPdeModifier<1>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 1.0,nullptr,0.0,1.0, 0.1));

            AbstractForce<1>* const p_force = new ChemoForceWithPdes<1>(0.1, 1e-4, 1e-4, 1e-4, p_pde_modifier);
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_force;
            delete p_force;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractForce<1>* p_force;
            input_arch >> p_force;

            TS_ASSERT_EQUALS(dynamic_cast<ChemoForceWithPdes<1>*>(p_force)->GetChemotacticSensitivity(), 0.1);

            delete p_force;
        }
    }

    void NoTestChemotacticForceWithPdeIn2D() 
    {
        // creation of the mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0u, false, 2, 5));
        nodes.push_back(new Node<2>(1u, false, 1, 5));
        nodes.push_back(new Node<2>(2u, false, 0, 5));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

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

        ChastePoint<2> lower(0.0, 0.0);
        ChastePoint<2> upper(10.0, 10.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        MAKE_PTR_ARGS(VegfEquationPde<2>, p_pde, (cell_population, 1.0, 1.0, 1.0, 0.1, 0.01));
        MAKE_PTR_ARGS(VegfBoundaryCondition<2>, p_bc, (1.0, 0));
        MAKE_PTR_ARGS(MolecularConcentrationsDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 1.0,nullptr,0.0,1.0, 0.1));
        p_pde_modifier->SetDependentVariableName("vegf_femesh_variable");
        p_pde_modifier->SetOutputGradient(true);
        p_pde_modifier->SetupInitialSolutionVector(cell_population);

        // we create the random force object 
        ChemoForceWithPdes<2> chemo_force_withpdes(0.1, 1e-4, 1e-4, 1e-4, p_pde_modifier);

        chemo_force_withpdes.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            if(p_cell->GetMutationState()->IsType<VesselTipMutationState>()){
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 1e-2, 1e-3);
            } else {
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.0, 1e-5); 
            }
        }

        // test archiving class 
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "chemo_force_withpdes.arch";
        {
            MAKE_PTR_ARGS(VegfEquationPde<2>, p_pde, (cell_population, 1.0, 1.0, 1.0, 0.1, 0.01));
            MAKE_PTR_ARGS(VegfBoundaryCondition<2>, p_bc, (1.0, 8));
            MAKE_PTR_ARGS(MolecularConcentrationsDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 1.0,nullptr,0.0,1.0, 0.1));

            AbstractForce<2>* const p_force = new ChemoForceWithPdes<2>(0.1, 1e-4, 1e-4, 1e-4, p_pde_modifier);
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

            TS_ASSERT_EQUALS(dynamic_cast<ChemoForceWithPdes<2>*>(p_force)->GetChemotacticSensitivity(), 0.1);

            delete p_force;
        }
    }

    void NoTestChemotacticForceWithPdeIn3D() 
    {
        // creation of the mesh
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u, false, 2, 5, 5));
        nodes.push_back(new Node<3>(1u, false, 1, 5, 5));
        nodes.push_back(new Node<3>(2u, false, 0, 5, 5));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        ChastePoint<3> lower(0.0, 0.0, 0.0);
        ChastePoint<3> upper(10.0, 10.0, 10.0);
        MAKE_PTR_ARGS(ChasteCuboid<3>, p_cuboid, (lower, upper));

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        MAKE_PTR_ARGS(VegfEquationPde<3>, p_pde, (cell_population, 1.0, 1.0, 1.0, 0.1, 0.01));
        MAKE_PTR_ARGS(VegfBoundaryCondition<3>, p_bc, (1.0, 0));
        MAKE_PTR_ARGS(MolecularConcentrationsDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 1.0,nullptr,0.0,1.0, 0.1));
        p_pde_modifier->SetDependentVariableName("vegf_femesh_variable");
        p_pde_modifier->SetOutputGradient(true);
        p_pde_modifier->SetupInitialSolutionVector(cell_population);

        // we create the random force object 
        ChemoForceWithPdes<3> chemo_force_withpdes(0.1, 1e-4, 1e-4, 1e-4, p_pde_modifier);

        chemo_force_withpdes.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            if(p_cell->GetMutationState()->IsType<VesselTipMutationState>()){
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 1e-2, 1e-3);
            } else {
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.0, 1e-5); 
            }
        }

        // test archiving class 
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "chemo_force_withpdes.arch";
        {
            MAKE_PTR_ARGS(VegfEquationPde<3>, p_pde, (cell_population, 1.0, 1.0, 1.0, 0.1, 0.01));
            MAKE_PTR_ARGS(VegfBoundaryCondition<3>, p_bc, (1.0, 8));
            MAKE_PTR_ARGS(MolecularConcentrationsDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 1.0,nullptr,0.0,1.0, 0.1));

            AbstractForce<3>* const p_force = new ChemoForceWithPdes<3>(0.1, 1e-4, 1e-4, 1e-4, p_pde_modifier);
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_force;
            delete p_force;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractForce<3>* p_force;
            input_arch >> p_force;

            TS_ASSERT_EQUALS(dynamic_cast<ChemoForceWithPdes<3>*>(p_force)->GetChemotacticSensitivity(), 0.1);

            delete p_force;
        }
    }

    void NoTestChemotacticForceWithApproximationOfPdeIn1D()
    {
        // creation of the mesh
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0u, false, 10));
        nodes.push_back(new Node<1>(1u, false, 9));
        nodes.push_back(new Node<1>(2u, false, 8));

        NodesOnlyMesh<1> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<1> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        // we create the random force object 
        ChemoForceWithAnalyticalApproximationPde<1> chemo_force_withanalyticalapproximationpde(1e-4, 1e-6, 1e4, 1.0, 0.0, 0.0, 0.5, 0.1);

        chemo_force_withanalyticalapproximationpde.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        { 
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            if(p_cell->GetMutationState()->IsType<VesselTipMutationState>()){
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[0], -3.7e-3, 1e-2); // value of gradient for the only tip cell in this configuration 
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 1e-5, 1e-5); 
            } else {
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 0.0);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.0, 1.0);
            }
        }

        // test archiving class 
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "chemo_force_withanalyticalapproximationpde.arch";
        {
            AbstractForce<1>* const p_force = new ChemoForceWithAnalyticalApproximationPde<1>(1e-4);
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_force;
            delete p_force;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractForce<1>* p_force;
            input_arch >> p_force;

            TS_ASSERT_EQUALS(dynamic_cast<ChemoForceWithAnalyticalApproximationPde<1>*>(p_force)->GetChemotacticSensitivity(), 0.1);

            delete p_force;
        }
    }

    void NoTestChemotacticForceWithApproximationOfPdeIn2D()
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
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

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
        ChemoForceWithAnalyticalApproximationPde<2> chemo_force_withanalyticalapproximationpde(1e-4, 1e-6, 1e4, 1.0, 0.0, 0.0, 0.5, 0.1);

        chemo_force_withanalyticalapproximationpde.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        { 
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            if(p_cell->GetMutationState()->IsType<VesselTipMutationState>()){
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[0], -3.7e-3, 1e-2); // value of gradient for the only tip cell in this configuration 
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 0.0);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 1e-5, 1e-5); // random norm
            } else {
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 0.0);
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 0.0);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.0, 1.0); // random norm
            }
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

            TS_ASSERT_EQUALS(dynamic_cast<ChemoForceWithAnalyticalApproximationPde<2>*>(p_force)->GetChemotacticSensitivity(), 0.1);

            delete p_force;
        }
    }

    void NoTestChemotacticForceWithApproximationOfPdeIn3D()
    {
        // creation of the mesh
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u, false, 10, 5, 5));
        nodes.push_back(new Node<3>(1u, false, 9, 5, 5));
        nodes.push_back(new Node<3>(2u, false, 8, 5, 5));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        // we create the random force object 
        ChemoForceWithAnalyticalApproximationPde<3> chemo_force_withanalyticalapproximationpde(1e-4, 1e-6, 1e4, 1.0, 0.0, 0.0, 0.5, 0.1);

        chemo_force_withanalyticalapproximationpde.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        { 
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            if(p_cell->GetMutationState()->IsType<VesselTipMutationState>()){
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[0], -3.7e-3, 1e-2); // value of gradient for the only tip cell in this configuration 
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 0.0);
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[2], 0.0);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 1e-5, 1e-5); 
            } else {
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 0.0);
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 0.0);
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[2], 0.0);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.0, 1.0); 
            }
        }

        // test archiving class 
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "chemo_force_withanalyticalapproximationpde.arch";
        {
            AbstractForce<3>* const p_force = new ChemoForceWithAnalyticalApproximationPde<3>(1e-4);
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_force;
            delete p_force;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractForce<3>* p_force;
            input_arch >> p_force;

            TS_ASSERT_EQUALS(dynamic_cast<ChemoForceWithAnalyticalApproximationPde<3>*>(p_force)->GetChemotacticSensitivity(), 0.1);

            delete p_force;
        }
    }

    void NoTestPersistenceForceIn1D()
    {
        // creation of the mesh
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0u, false, 10));
        nodes.push_back(new Node<1>(1u, false, 9));
        nodes.push_back(new Node<1>(2u, false, 8));

        NodesOnlyMesh<1> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<1> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        cell_population.Update();

        std::vector<std::string> mDataLabels = std::vector<std::string>{ "fx", "fy", "fz" };
        for ( typename NodeBasedCellPopulation<1>::Iterator cell_it = cell_population.Begin(); 
                                                cell_it != cell_population.End(); ++cell_it )
        {
            for (unsigned i = 0; i < 1; i++)
            {
                (*cell_it)->GetCellData()->SetItem(mDataLabels[i], 0.0);
            }
        }

        // we create the random force object 
        PersistenceForce<1> persistence_force(0.1);
        persistence_force.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            if(p_cell->GetMutationState()->IsType<VesselTipMutationState>()){
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 0.1, 1e-1);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.1, 1e-1);  
            } else {
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 0.0);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.0, 1.0); 
            }
        }

        // test archiving class 
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "persistence_force.arch";
        {
            AbstractForce<1>* const p_force = new PersistenceForce<1>(0.1);
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_force;
            delete p_force;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractForce<1>* p_force;
            input_arch >> p_force;

            TS_ASSERT_EQUALS(dynamic_cast<PersistenceForce<1>*>(p_force)->GetPersistenceCoefficient(), 0.1);

            delete p_force;
        }
    }

    void NoTestPersistenceForceIn2D()
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
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

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

        cell_population.Update();

        std::vector<std::string> mDataLabels = std::vector<std::string>{ "fx", "fy", "fz" };
        for ( typename NodeBasedCellPopulation<2>::Iterator cell_it = cell_population.Begin(); 
                                                cell_it != cell_population.End(); ++cell_it )
        {
            for (unsigned i = 0; i < 2; i++)
            {
                (*cell_it)->GetCellData()->SetItem(mDataLabels[i], 0.0);
            }
        }

        // we create the random force object 
        PersistenceForce<2> persistence_force(0.1);
        persistence_force.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            if(p_cell->GetMutationState()->IsType<VesselTipMutationState>()){
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 0.1, 1e-1);
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 0.1, 1e-1);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.1, 1e-1);  
            } else {
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 0.0);
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 0.0);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.0, 1.0); 
            }
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

    void NoTestPersistenceForceIn3D()
    {
        // creation of the mesh
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u, false, 10, 5, 5));
        nodes.push_back(new Node<3>(1u, false, 9, 5, 5));
        nodes.push_back(new Node<3>(2u, false, 8, 5, 5));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        cell_population.Update();

        std::vector<std::string> mDataLabels = std::vector<std::string>{ "fx", "fy", "fz" };
        for ( typename NodeBasedCellPopulation<3>::Iterator cell_it = cell_population.Begin(); 
                                                cell_it != cell_population.End(); ++cell_it )
        {
            for (unsigned i = 0; i < 3; i++)
            {
                (*cell_it)->GetCellData()->SetItem(mDataLabels[i], 0.0);
            }
        }

        // we create the random force object 
        PersistenceForce<3> persistence_force(0.1);
        persistence_force.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            if(p_cell->GetMutationState()->IsType<VesselTipMutationState>()){
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 0.1, 1e-1);
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 0.1, 1e-1);
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[2], 0.1, 1e-1);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.1, 1e-1);  
            } else {
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 0.0);
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 0.0);
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[2], 0.0);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.0, 1.0); 
            }
        }

        // test archiving class 
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "persistence_force.arch";
        {
            AbstractForce<3>* const p_force = new PersistenceForce<3>(0.1);
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_force;
            delete p_force;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractForce<3>* p_force;
            input_arch >> p_force;

            TS_ASSERT_EQUALS(dynamic_cast<PersistenceForce<3>*>(p_force)->GetPersistenceCoefficient(), 0.1);

            delete p_force;
        }
    }

    void NoTestAngularForceIn1D()
    {
        // creation of the mesh
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0u, false, 10));
        nodes.push_back(new Node<1>(1u, false, 9));
        nodes.push_back(new Node<1>(2u, false, 8));

        NodesOnlyMesh<1> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states 
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<1> cell_population(mesh, cells);

        // fully constrain the first cell using the boundary condition 
        unsigned node_index_tip_cell = cell_population.GetLocationIndexUsingCell(0);
        std::vector<unsigned> pinned_node_indices;
        pinned_node_indices.push_back(node_index_tip_cell);
        typedef PinnedCellsBoundaryCondition<1,1> PinnedCellsBoundaryCondition;
        MAKE_PTR_ARGS(PinnedCellsBoundaryCondition, p_boundary_condition, (&cell_population, pinned_node_indices));

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        cell_population.Update();

        // Set up cell data on the cell population : initialisation for the division 
        for (typename NodeBasedCellPopulation<1>::Iterator cell_it = cell_population.Begin(); 
                                                cell_it != cell_population.End(); ++cell_it )
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_it);
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            c_vector<double, 1> node_coordinates = cell_population.GetLocationOfCellCentre(p_cell);

            double tortuosity = norm_2(node_coordinates);

            double t = SimulationTime::Instance()->GetTime(); 

            (*cell_it)->GetCellData()->SetItem("node_index", node_index);
            (*cell_it)->GetCellData()->SetItem("daughter_type", 0.0); // 0.0 = no daughter cell, 1.0 = daughter cell
            (*cell_it)->GetCellData()->SetItem("BranchNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingPoint", 0.0); // 0.0 = no, 1.0 = yes : necessary to modify the mutation state after sprouting in this modifier 
            (*cell_it)->GetCellData()->SetItem("DivisionNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingLeader", 0.0); // indice of the branching cell leading the branch of the cell considered 
            (*cell_it)->GetCellData()->SetItem("BranchingSegment", 0.0);
            (*cell_it)->GetCellData()->SetItem("ConsecutiveBranches", 1.0); 
            (*cell_it)->GetCellData()->SetItem("BirthTimeCell", t);
            (*cell_it)->GetCellData()->SetItem("TortuosityParent", tortuosity);
            (*cell_it)->GetCellData()->SetItem("TortuosityDaughter", tortuosity);
            (*cell_it)->GetCellData()->SetItem("LoopNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("Anastomosis", 0);
        }

        // we create the random force object 
        AngularForce<1> angular_force(1e-5);
        angular_force.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            if(p_cell->GetMutationState()->IsType<VesselSegmentMutationState>()){
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 1e-5, 1e-5);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 1e-5, 1e-5);  
            } else {
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 0.0);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.0, 1.0); 
            }
        }

        // test archiving class 
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "angular_force.arch";
        {
            AbstractForce<1>* const p_force = new AngularForce<1>(1e-5);
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_force;
            delete p_force;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractForce<1>* p_force;
            input_arch >> p_force;

            TS_ASSERT_EQUALS(dynamic_cast<AngularForce<1>*>(p_force)->GetAngularPersistence(), 1e-5);

            delete p_force;
        }
    }

    void NoTestAngularForceIn2D()
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
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

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

        // fully constrain the first cell using the boundary condition 
        unsigned node_index_tip_cell = cell_population.GetLocationIndexUsingCell(0);
        std::vector<unsigned> pinned_node_indices;
        pinned_node_indices.push_back(node_index_tip_cell);
        typedef PinnedCellsBoundaryCondition<2,2> PinnedCellsBoundaryCondition;
        MAKE_PTR_ARGS(PinnedCellsBoundaryCondition, p_boundary_condition, (&cell_population, pinned_node_indices));

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        cell_population.Update();

        // Set up cell data on the cell population : initialisation for the division 
        for (typename NodeBasedCellPopulation<2>::Iterator cell_it = cell_population.Begin(); 
                                                cell_it != cell_population.End(); ++cell_it )
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_it);
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            c_vector<double, 2> node_coordinates = cell_population.GetLocationOfCellCentre(p_cell);

            double tortuosity = norm_2(node_coordinates);

            double t = SimulationTime::Instance()->GetTime(); 

            (*cell_it)->GetCellData()->SetItem("node_index", node_index);
            (*cell_it)->GetCellData()->SetItem("daughter_type", 0.0); // 0.0 = no daughter cell, 1.0 = daughter cell
            (*cell_it)->GetCellData()->SetItem("BranchNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingPoint", 0.0); // 0.0 = no, 1.0 = yes : necessary to modify the mutation state after sprouting in this modifier 
            (*cell_it)->GetCellData()->SetItem("DivisionNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingLeader", 0.0); // indice of the branching cell leading the branch of the cell considered 
            (*cell_it)->GetCellData()->SetItem("BranchingSegment", 0.0);
            (*cell_it)->GetCellData()->SetItem("ConsecutiveBranches", 1.0); 
            (*cell_it)->GetCellData()->SetItem("BirthTimeCell", t);
            (*cell_it)->GetCellData()->SetItem("TortuosityParent", tortuosity);
            (*cell_it)->GetCellData()->SetItem("TortuosityDaughter", tortuosity);
            (*cell_it)->GetCellData()->SetItem("LoopNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("Anastomosis", 0);
        }

        // we create the random force object 
        AngularForce<2> angular_force(1e-5);

        angular_force.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            if(p_cell->GetMutationState()->IsType<VesselSegmentMutationState>()){
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 1e-5, 1e-5);
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 1e-5, 1e-5);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 1e-5, 1e-5);  
            } else {
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 0.0);
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 0.0);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.0, 1.0); 
            }
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

    void NoTestAngularForceIn3D()
    {
        // creation of the mesh
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u, false, 10, 5, 5));
        nodes.push_back(new Node<3>(1u, false, 9, 5, 5));
        nodes.push_back(new Node<3>(2u, false, 8, 5, 5));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states 
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        // fully constrain the first cell using the boundary condition 
        unsigned node_index_tip_cell = cell_population.GetLocationIndexUsingCell(0);
        std::vector<unsigned> pinned_node_indices;
        pinned_node_indices.push_back(node_index_tip_cell);
        typedef PinnedCellsBoundaryCondition<3,3> PinnedCellsBoundaryCondition;
        MAKE_PTR_ARGS(PinnedCellsBoundaryCondition, p_boundary_condition, (&cell_population, pinned_node_indices));

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        cell_population.Update();

        // Set up cell data on the cell population : initialisation for the division 
        for (typename NodeBasedCellPopulation<3>::Iterator cell_it = cell_population.Begin(); 
                                                cell_it != cell_population.End(); ++cell_it )
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_it);
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            c_vector<double, 3> node_coordinates = cell_population.GetLocationOfCellCentre(p_cell);

            double tortuosity = norm_2(node_coordinates);

            double t = SimulationTime::Instance()->GetTime(); 

            (*cell_it)->GetCellData()->SetItem("node_index", node_index);
            (*cell_it)->GetCellData()->SetItem("daughter_type", 0.0); // 0.0 = no daughter cell, 1.0 = daughter cell
            (*cell_it)->GetCellData()->SetItem("BranchNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingPoint", 0.0); // 0.0 = no, 1.0 = yes : necessary to modify the mutation state after sprouting in this modifier 
            (*cell_it)->GetCellData()->SetItem("DivisionNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingLeader", 0.0); // indice of the branching cell leading the branch of the cell considered 
            (*cell_it)->GetCellData()->SetItem("BranchingSegment", 0.0);
            (*cell_it)->GetCellData()->SetItem("ConsecutiveBranches", 1.0); 
            (*cell_it)->GetCellData()->SetItem("BirthTimeCell", t);
            (*cell_it)->GetCellData()->SetItem("TortuosityParent", tortuosity);
            (*cell_it)->GetCellData()->SetItem("TortuosityDaughter", tortuosity);
            (*cell_it)->GetCellData()->SetItem("LoopNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("Anastomosis", 0);
        }

        // we create the random force object 
        AngularForce<3> angular_force(1e-5);

        angular_force.AddForceContribution(cell_population);

        // test object 
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            if(p_cell->GetMutationState()->IsType<VesselSegmentMutationState>()){
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 1e-5, 1e-5);
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 1e-5, 1e-5);
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[2], 1e-5, 1e-5);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 1e-5, 1e-5);  
            } else {
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 0.0);
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 0.0);
                TS_ASSERT_EQUALS(cell_population.GetNode(node_index)->rGetAppliedForce()[2], 0.0);
                TS_ASSERT_DELTA(norm_2(cell_population.GetNode(node_index)->rGetAppliedForce()), 0.0, 1.0); 
            }
        }

        // test archiving class 
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "angular_force.arch";
        {
            AbstractForce<3>* const p_force = new AngularForce<3>(1e-5);
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_force;
            delete p_force;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractForce<3>* p_force;
            input_arch >> p_force;

            TS_ASSERT_EQUALS(dynamic_cast<AngularForce<3>*>(p_force)->GetAngularPersistence(), 1e-5);

            delete p_force;
        }
    }

    void NoTestPinnedCellBoundaryConditionIn1D()
    {
        // creation of the mesh
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0u, false, 10));
        nodes.push_back(new Node<1>(1u, false, 9));
        nodes.push_back(new Node<1>(2u, false, 8));

        NodesOnlyMesh<1> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<1> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        // create the cell population boundary condition object 
        std::vector<unsigned> pinned_node_indices;
        // fully constrain the first cell using the boundary condition 
        unsigned node_index_tip_cell = cell_population.GetLocationIndexUsingCell(0);
        pinned_node_indices.push_back(node_index_tip_cell);
        PinnedCellsBoundaryCondition<1,1> bc(&cell_population, pinned_node_indices);

        // we store the previous node locations 
        std::map<Node<1>*, c_vector<double, 1> > old_node_locations;
        for (AbstractMesh<1,1>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                node_iter != mesh.GetNodeIteratorEnd();
                ++node_iter)
        {
            old_node_locations[&(*node_iter)] = node_iter->rGetLocation();
        }
        bc.ImposeBoundaryCondition(old_node_locations);

        // now we check that the cell population verifies the boundary condition
        bool population_satisfies_bc = bc.VerifyBoundaryCondition();
        TS_ASSERT_EQUALS(population_satisfies_bc, true);

        // move node and check return false 
        c_vector<double,1> new_node_position;
        new_node_position(0) = 9;
        ChastePoint<1> new_node(new_node_position);
        cell_population.SetNode(node_index_tip_cell, new_node);
        cell_population.Update();
        population_satisfies_bc = bc.VerifyBoundaryCondition();
        TS_ASSERT_EQUALS(population_satisfies_bc, false);
    }

    void NoTestPinnedCellBoundaryConditionIn2D()
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
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

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

        // create the cell population boundary condition object 
        std::vector<unsigned> pinned_node_indices;
        // fully constrain the first cell using the boundary condition 
        unsigned node_index_tip_cell = cell_population.GetLocationIndexUsingCell(0);
        pinned_node_indices.push_back(node_index_tip_cell);
        PinnedCellsBoundaryCondition<2,2> bc(&cell_population, pinned_node_indices);

        // we store the previous node locations 
        std::map<Node<2>*, c_vector<double, 2> > old_node_locations;
        for (AbstractMesh<2,2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                node_iter != mesh.GetNodeIteratorEnd();
                ++node_iter)
        {
            old_node_locations[&(*node_iter)] = node_iter->rGetLocation();
        }
        bc.ImposeBoundaryCondition(old_node_locations);

        // now we check that the cell population verifies the boundary condition
        bool population_satisfies_bc = bc.VerifyBoundaryCondition();
        TS_ASSERT_EQUALS(population_satisfies_bc, true);

        // move node and check return false 
        c_vector<double,2> new_node_position;
        new_node_position(0) = 9;
        new_node_position(1) = 4;
        ChastePoint<2> new_node(new_node_position);
        cell_population.SetNode(node_index_tip_cell, new_node);
        cell_population.Update();
        population_satisfies_bc = bc.VerifyBoundaryCondition();
        TS_ASSERT_EQUALS(population_satisfies_bc, false);
    }

    void NoTestPinnedCellBoundaryConditionIn3D()
    {
        // creation of the mesh
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u, false, 10, 5, 5));
        nodes.push_back(new Node<3>(1u, false, 9, 5, 5));
        nodes.push_back(new Node<3>(2u, false, 8, 5, 5));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        // create the cell population boundary condition object 
        std::vector<unsigned> pinned_node_indices;
        // fully constrain the first cell using the boundary condition 
        unsigned node_index_tip_cell = cell_population.GetLocationIndexUsingCell(0);
        pinned_node_indices.push_back(node_index_tip_cell);
        PinnedCellsBoundaryCondition<3,3> bc(&cell_population, pinned_node_indices);

        // we store the previous node locations 
        std::map<Node<3>*, c_vector<double, 3> > old_node_locations;
        for (AbstractMesh<3,3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                node_iter != mesh.GetNodeIteratorEnd();
                ++node_iter)
        {
            old_node_locations[&(*node_iter)] = node_iter->rGetLocation();
        }
        bc.ImposeBoundaryCondition(old_node_locations);

        // now we check that the cell population verifies the boundary condition
        bool population_satisfies_bc = bc.VerifyBoundaryCondition();
        TS_ASSERT_EQUALS(population_satisfies_bc, true);

        // move node and check return false 
        c_vector<double,3> new_node_position;
        new_node_position(0) = 9;
        new_node_position(1) = 4;
        new_node_position(2) = 4;
        ChastePoint<3> new_node(new_node_position);
        cell_population.SetNode(node_index_tip_cell, new_node);
        cell_population.Update();
        population_satisfies_bc = bc.VerifyBoundaryCondition();
        TS_ASSERT_EQUALS(population_satisfies_bc, false);
    }

    void NoTestSproutingRuleIn1D() 
    {
        // creation of the mesh
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0u, false, 10));
        nodes.push_back(new Node<1>(1u, false, 9));
        nodes.push_back(new Node<1>(2u, false, 8));

        NodesOnlyMesh<1> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<1> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        cell_population.Update();

        // Set up cell data on the cell population : initialisation for the division 
        for (typename NodeBasedCellPopulation<1>::Iterator cell_it = cell_population.Begin(); 
                                                cell_it != cell_population.End(); ++cell_it )
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_it);
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            c_vector<double, 1> node_coordinates = cell_population.GetLocationOfCellCentre(p_cell);

            double tortuosity = norm_2(node_coordinates);

            double t = SimulationTime::Instance()->GetTime(); 

            (*cell_it)->GetCellData()->SetItem("node_index", node_index);
            (*cell_it)->GetCellData()->SetItem("daughter_type", 0.0); // 0.0 = no daughter cell, 1.0 = daughter cell
            (*cell_it)->GetCellData()->SetItem("BranchNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingPoint", 0.0); // 0.0 = no, 1.0 = yes : necessary to modify the mutation state after sprouting in this modifier 
            (*cell_it)->GetCellData()->SetItem("DivisionNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingLeader", 0.0); // indice of the branching cell leading the branch of the cell considered 
            (*cell_it)->GetCellData()->SetItem("BranchingSegment", 0.0);
            (*cell_it)->GetCellData()->SetItem("ConsecutiveBranches", 1.0); 
            (*cell_it)->GetCellData()->SetItem("BirthTimeCell", t);
            (*cell_it)->GetCellData()->SetItem("TortuosityParent", tortuosity);
            (*cell_it)->GetCellData()->SetItem("TortuosityDaughter", tortuosity);
            (*cell_it)->GetCellData()->SetItem("LoopNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("Anastomosis", 0);
        }

        // create a vector of the expected coordinates of the parent and daughter cells : in this case, the tip cell divides into a vessel segment
        CellPtr p_cell0 = cell_population.GetCellUsingLocationIndex(2);
        c_vector<double,1> old_parent_location = cell_population.GetLocationOfCellCentre(p_cell0);

        // 0.3 is the MeinekeDivisionSeparation value and (-1, 0) is the daughter_direction
        c_vector<double, 1> expected_daughter_location_vessel_segment  = old_parent_location;
        c_vector<double, 1> expected_parent_location_vessel_segment;
        expected_parent_location_vessel_segment[0] = old_parent_location[0] - 0.3*1;  

        c_vector<double, 2> expected_daughter_location_tip_cell;
        expected_daughter_location_tip_cell[0] = old_parent_location[0] - 0.5*0.3*(-1) + 0.5*0.3*1; 
        c_vector<double, 2> expected_parent_location_tip_cell;
        expected_parent_location_tip_cell[0] = old_parent_location[0] + 0.5*0.3*(-1) + 0.5*0.3*1; 

        // we create the sprouting rule object with a sprouting proba of 0 (to obtain a division into a vessel segment) but the parent class does not have rGetDaughterLocation implemented 
        typedef SproutingRule<1,1> SproutingRule;
        MAKE_PTR_ARGS(SproutingRule, p_sprouting_rule_vessel_segment, (0.01, 2.0));
        // TS_ASSERT_DELTA(p_sprouting_rule_vessel_segment->rGetDaughterLocation()[0], 7.7, 1e-1);

        MAKE_PTR_ARGS(SproutingRule, p_sprouting_rule_tip_cell, (0.99, 2.0));
        // TS_ASSERT_DELTA(p_sprouting_rule_tip_cell->rGetDaughterLocation()[0], 8.15, 1e-2);

        // case of vessel segment 
        // set the division rule to be the division rule for the cell population
        cell_population.SetCentreBasedDivisionRule(p_sprouting_rule_vessel_segment);
        boost::shared_ptr<AbstractCentreBasedDivisionRule<1,1> > p_sprouting_rule_vessel_segment_test = cell_population.GetCentreBasedDivisionRule();

        // test to check that the division rule returns the correct pair of vectors 
        std::pair<c_vector<double, 1>, c_vector<double, 1> > positions_vessel_segment = p_sprouting_rule_vessel_segment_test->CalculateCellDivisionVector(p_cell0, cell_population);

        c_vector<double, 1> parent_location_vessel_segment;
        parent_location_vessel_segment = positions_vessel_segment.first;
        TS_ASSERT_DELTA(parent_location_vessel_segment[0], expected_parent_location_vessel_segment[0], 1e-6);

        c_vector<double, 1> daughter_location_vessel_segment;
        daughter_location_vessel_segment = positions_vessel_segment.second;
        TS_ASSERT_DELTA(daughter_location_vessel_segment[0], expected_daughter_location_vessel_segment[0], 1e-6);

        // case of tip cell (take into account new position of the tip cell) 
        // set the division rule to be the division rule for the cell population
        cell_population.SetCentreBasedDivisionRule(p_sprouting_rule_tip_cell);
        boost::shared_ptr<AbstractCentreBasedDivisionRule<1,1> > p_sprouting_rule_tip_cell_test = cell_population.GetCentreBasedDivisionRule();

        // test to check that the division rule returns the correct pair of vectors 
        std::pair<c_vector<double, 1>, c_vector<double, 1> > positions_tip_cell = p_sprouting_rule_tip_cell_test->CalculateCellDivisionVector(p_cell0, cell_population);

        c_vector<double, 1> parent_location_tip_cell;
        parent_location_tip_cell = positions_tip_cell.first;
        TS_ASSERT_DELTA(parent_location_tip_cell[0], expected_parent_location_tip_cell[0], 1e-6);

        c_vector<double, 1> daughter_location_tip_cell;
        daughter_location_tip_cell = positions_tip_cell.second;
        TS_ASSERT_DELTA(daughter_location_tip_cell[0], expected_daughter_location_tip_cell[0], 1e-6);

        // archiving
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "sprouting_rule.arch";
        {
            SproutingRule* const p_sprouting_rule = new SproutingRule(0.08, 2.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_sprouting_rule;

            delete p_sprouting_rule;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            SproutingRule* p_sprouting_rule;
            input_arch >> p_sprouting_rule;

            TS_ASSERT_EQUALS(dynamic_cast<SproutingRule*>(p_sprouting_rule)->GetMaxSproutingRate(), 0.08);

            delete p_sprouting_rule;
        }
    }

    void NoTestSproutingRuleIn2D() 
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
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

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

        cell_population.Update();

        // Set up cell data on the cell population : initialisation for the division 
        for (typename NodeBasedCellPopulation<2>::Iterator cell_it = cell_population.Begin(); 
                                                cell_it != cell_population.End(); ++cell_it )
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_it);
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            c_vector<double, 2> node_coordinates = cell_population.GetLocationOfCellCentre(p_cell);

            double tortuosity = norm_2(node_coordinates);

            double t = SimulationTime::Instance()->GetTime(); 

            (*cell_it)->GetCellData()->SetItem("node_index", node_index);
            (*cell_it)->GetCellData()->SetItem("daughter_type", 0.0); // 0.0 = no daughter cell, 1.0 = daughter cell
            (*cell_it)->GetCellData()->SetItem("BranchNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingPoint", 0.0); // 0.0 = no, 1.0 = yes : necessary to modify the mutation state after sprouting in this modifier 
            (*cell_it)->GetCellData()->SetItem("DivisionNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingLeader", 0.0); // indice of the branching cell leading the branch of the cell considered 
            (*cell_it)->GetCellData()->SetItem("BranchingSegment", 0.0);
            (*cell_it)->GetCellData()->SetItem("ConsecutiveBranches", 1.0); 
            (*cell_it)->GetCellData()->SetItem("BirthTimeCell", t);
            (*cell_it)->GetCellData()->SetItem("TortuosityParent", tortuosity);
            (*cell_it)->GetCellData()->SetItem("TortuosityDaughter", tortuosity);
            (*cell_it)->GetCellData()->SetItem("LoopNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("Anastomosis", 0);
        }

        // create a vector of the expected coordinates of the parent and daughter cells : in this case, the tip cell divides into a vessel segment
        CellPtr p_cell0 = cell_population.GetCellUsingLocationIndex(2);
        c_vector<double,2> old_parent_location = cell_population.GetLocationOfCellCentre(p_cell0);

        // 0.3 is the MeinekeDivisionSeparation value and (-1, 0) is the daughter_direction
        c_vector<double, 2> expected_daughter_location_vessel_segment  = old_parent_location;
        c_vector<double, 2> expected_parent_location_vessel_segment;
        expected_parent_location_vessel_segment[0] = old_parent_location[0] - 0.3*1;  
        expected_parent_location_vessel_segment[1] = old_parent_location[1] - 0.3*0;

        c_vector<double, 2> expected_daughter_location_tip_cell;
        expected_daughter_location_tip_cell[0] = old_parent_location[0] - 0.5*0.3*0 + 0.5*0.3*1;  
        expected_daughter_location_tip_cell[1] = old_parent_location[1] - 0.5*0.3*1 + 0.5*0.3*0;
        c_vector<double, 2> expected_parent_location_tip_cell;
        expected_parent_location_tip_cell[0] = old_parent_location[0] + 0.5*0.3*0 + 0.5*0.3*1;  
        expected_parent_location_tip_cell[1] = old_parent_location[1] + 0.5*0.3*1 + 0.5*0.3*0;

        // we create the sprouting rule object with a sprouting proba of 0 (to obtain a division into a vessel segment)
        typedef SproutingRule<2,2> SproutingRule;
        MAKE_PTR_ARGS(SproutingRule, p_sprouting_rule_vessel_segment, (0.01, 2.0));
        // TS_ASSERT_DELTA(p_sprouting_rule_vessel_segment->rGetDaughterLocation()[0], 7.7, 1e-1);
        // TS_ASSERT_DELTA(p_sprouting_rule_vessel_segment->rGetDaughterLocation()[1], 0.0, 1e-1);

        MAKE_PTR_ARGS(SproutingRule, p_sprouting_rule_tip_cell, (0.99, 2.0));
        // TS_ASSERT_DELTA(p_sprouting_rule_tip_cell->rGetDaughterLocation()[0], 8.15, 1e-2);
        // TS_ASSERT_DELTA(p_sprouting_rule_tip_cell->rGetDaughterLocation()[1], 0.15, 1e-2);

        // case of vessel segment 
        // set the division rule to be the division rule for the cell population
        cell_population.SetCentreBasedDivisionRule(p_sprouting_rule_vessel_segment);
        boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_sprouting_rule_vessel_segment_test = cell_population.GetCentreBasedDivisionRule();

        // test to check that the division rule returns the correct pair of vectors 
        std::pair<c_vector<double, 2>, c_vector<double, 2> > positions_vessel_segment = p_sprouting_rule_vessel_segment_test->CalculateCellDivisionVector(p_cell0, cell_population);

        c_vector<double, 2> parent_location_vessel_segment;
        parent_location_vessel_segment = positions_vessel_segment.first;
        TS_ASSERT_DELTA(parent_location_vessel_segment[0], expected_parent_location_vessel_segment[0], 1e-6);
        TS_ASSERT_DELTA(parent_location_vessel_segment[1], expected_parent_location_vessel_segment[1], 1e-6);

        c_vector<double, 2> daughter_location_vessel_segment;
        daughter_location_vessel_segment = positions_vessel_segment.second;
        TS_ASSERT_DELTA(daughter_location_vessel_segment[0], expected_daughter_location_vessel_segment[0], 1e-6);
        TS_ASSERT_DELTA(daughter_location_vessel_segment[1], expected_daughter_location_vessel_segment[1], 1e-6);

        // case of tip cell (take into account new position of the tip cell) 
        // set the division rule to be the division rule for the cell population
        cell_population.SetCentreBasedDivisionRule(p_sprouting_rule_tip_cell);
        boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_sprouting_rule_tip_cell_test = cell_population.GetCentreBasedDivisionRule();

        // test to check that the division rule returns the correct pair of vectors 
        std::pair<c_vector<double, 2>, c_vector<double, 2> > positions_tip_cell = p_sprouting_rule_tip_cell_test->CalculateCellDivisionVector(p_cell0, cell_population);

        c_vector<double, 2> parent_location_tip_cell;
        parent_location_tip_cell = positions_tip_cell.first;
        TS_ASSERT_DELTA(parent_location_tip_cell[0], expected_parent_location_tip_cell[0], 1e-6);
        TS_ASSERT_DELTA(parent_location_tip_cell[1], expected_parent_location_tip_cell[1], 1e-6);

        c_vector<double, 2> daughter_location_tip_cell;
        daughter_location_tip_cell = positions_tip_cell.second;
        TS_ASSERT_DELTA(daughter_location_tip_cell[0], expected_daughter_location_tip_cell[0], 1e-6);
        TS_ASSERT_DELTA(daughter_location_tip_cell[1], expected_daughter_location_tip_cell[1], 1e-6);

        // archiving
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "sprouting_rule.arch";
        {
            SproutingRule* const p_sprouting_rule = new SproutingRule(0.08, 2.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_sprouting_rule;

            delete p_sprouting_rule;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            SproutingRule* p_sprouting_rule;
            input_arch >> p_sprouting_rule;

            TS_ASSERT_EQUALS(dynamic_cast<SproutingRule*>(p_sprouting_rule)->GetMaxSproutingRate(), 0.08);

            delete p_sprouting_rule;
        }
    }

    void NoTestSproutingRuleIn3D() 
    {
        // creation of the mesh
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u, false, 10, 5, 5));
        nodes.push_back(new Node<3>(1u, false, 9, 5, 5));
        nodes.push_back(new Node<3>(2u, false, 8, 5, 5));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        cell_population.Update();

        // Set up cell data on the cell population : initialisation for the division 
        for (typename NodeBasedCellPopulation<3>::Iterator cell_it = cell_population.Begin(); 
                                                cell_it != cell_population.End(); ++cell_it )
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_it);
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            c_vector<double, 3> node_coordinates = cell_population.GetLocationOfCellCentre(p_cell);

            double tortuosity = norm_2(node_coordinates);

            double t = SimulationTime::Instance()->GetTime(); 

            (*cell_it)->GetCellData()->SetItem("node_index", node_index);
            (*cell_it)->GetCellData()->SetItem("daughter_type", 0.0); // 0.0 = no daughter cell, 1.0 = daughter cell
            (*cell_it)->GetCellData()->SetItem("BranchNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingPoint", 0.0); // 0.0 = no, 1.0 = yes : necessary to modify the mutation state after sprouting in this modifier 
            (*cell_it)->GetCellData()->SetItem("DivisionNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingLeader", 0.0); // indice of the branching cell leading the branch of the cell considered 
            (*cell_it)->GetCellData()->SetItem("BranchingSegment", 0.0);
            (*cell_it)->GetCellData()->SetItem("ConsecutiveBranches", 1.0); 
            (*cell_it)->GetCellData()->SetItem("BirthTimeCell", t);
            (*cell_it)->GetCellData()->SetItem("TortuosityParent", tortuosity);
            (*cell_it)->GetCellData()->SetItem("TortuosityDaughter", tortuosity);
            (*cell_it)->GetCellData()->SetItem("LoopNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("Anastomosis", 0);
        }

        // create a vector of the expected coordinates of the parent and daughter cells : in this case, the tip cell divides into a vessel segment
        CellPtr p_cell0 = cell_population.GetCellUsingLocationIndex(2);
        c_vector<double,3> old_parent_location = cell_population.GetLocationOfCellCentre(p_cell0);

        // 0.3 is the MeinekeDivisionSeparation value and (-1, 0) is the daughter_direction
        c_vector<double, 3> expected_daughter_location_vessel_segment  = old_parent_location;
        c_vector<double, 3> expected_parent_location_vessel_segment;
        expected_parent_location_vessel_segment[0] = old_parent_location[0] - 0.3*1;  
        expected_parent_location_vessel_segment[1] = old_parent_location[1] - 0.3*0;
        expected_parent_location_vessel_segment[2] = old_parent_location[2] - 0.3*0;

        c_vector<double, 3> expected_daughter_location_tip_cell;
        expected_daughter_location_tip_cell[0] = old_parent_location[0] - 0.5*0.3*0 + 0.5*0.3*1;  
        expected_daughter_location_tip_cell[1] = old_parent_location[1] - 0.5*0.3*1 + 0.5*0.3*0;
        expected_daughter_location_tip_cell[2] = old_parent_location[2] - 0.5*0.3*0 + 0.5*0.3*0;
        c_vector<double, 3> expected_parent_location_tip_cell;
        expected_parent_location_tip_cell[0] = old_parent_location[0] + 0.5*0.3*0 + 0.5*0.3*1;  
        expected_parent_location_tip_cell[1] = old_parent_location[1] + 0.5*0.3*1 + 0.5*0.3*0;
        expected_parent_location_tip_cell[2] = old_parent_location[2] + 0.5*0.3*0 + 0.5*0.3*0;

        // we create the sprouting rule object with a sprouting proba of 0 (to obtain a division into a vessel segment)
        typedef SproutingRule<3,3> SproutingRule;
        MAKE_PTR_ARGS(SproutingRule, p_sprouting_rule_vessel_segment, (0.01, 2.0));
        // TS_ASSERT_DELTA(p_sprouting_rule_vessel_segment->rGetDaughterLocation()[0], 7.7, 1e-1);
        // TS_ASSERT_DELTA(p_sprouting_rule_vessel_segment->rGetDaughterLocation()[1], 0.0, 1e-1);
        // TS_ASSERT_DELTA(p_sprouting_rule_vessel_segment->rGetDaughterLocation()[2], 0.0, 1e-1);

        MAKE_PTR_ARGS(SproutingRule, p_sprouting_rule_tip_cell, (0.99, 2.0));
        // TS_ASSERT_DELTA(p_sprouting_rule_tip_cell->rGetDaughterLocation()[0], 8.15, 1e-2);
        // TS_ASSERT_DELTA(p_sprouting_rule_tip_cell->rGetDaughterLocation()[1], 0.15, 1e-2);
        // TS_ASSERT_DELTA(p_sprouting_rule_tip_cell->rGetDaughterLocation()[2], 0.15, 1e-2);

        // case of vessel segment 
        // set the division rule to be the division rule for the cell population
        cell_population.SetCentreBasedDivisionRule(p_sprouting_rule_vessel_segment);
        boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_sprouting_rule_vessel_segment_test = cell_population.GetCentreBasedDivisionRule();

        // test to check that the division rule returns the correct pair of vectors 
        std::pair<c_vector<double, 3>, c_vector<double, 3> > positions_vessel_segment = p_sprouting_rule_vessel_segment_test->CalculateCellDivisionVector(p_cell0, cell_population);

        c_vector<double, 3> parent_location_vessel_segment;
        parent_location_vessel_segment = positions_vessel_segment.first;
        TS_ASSERT_DELTA(parent_location_vessel_segment[0], expected_parent_location_vessel_segment[0], 1e-6);
        TS_ASSERT_DELTA(parent_location_vessel_segment[1], expected_parent_location_vessel_segment[1], 1e-6);
        TS_ASSERT_DELTA(parent_location_vessel_segment[2], expected_parent_location_vessel_segment[2], 1e-6);

        c_vector<double, 3> daughter_location_vessel_segment;
        daughter_location_vessel_segment = positions_vessel_segment.second;
        TS_ASSERT_DELTA(daughter_location_vessel_segment[0], expected_daughter_location_vessel_segment[0], 1e-6);
        TS_ASSERT_DELTA(daughter_location_vessel_segment[1], expected_daughter_location_vessel_segment[1], 1e-6);
        TS_ASSERT_DELTA(daughter_location_vessel_segment[2], expected_daughter_location_vessel_segment[2], 1e-6);

        // case of tip cell (take into account new position of the tip cell) 
        // set the division rule to be the division rule for the cell population
        cell_population.SetCentreBasedDivisionRule(p_sprouting_rule_tip_cell);
        boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_sprouting_rule_tip_cell_test = cell_population.GetCentreBasedDivisionRule();

        // test to check that the division rule returns the correct pair of vectors 
        std::pair<c_vector<double, 3>, c_vector<double, 3> > positions_tip_cell = p_sprouting_rule_tip_cell_test->CalculateCellDivisionVector(p_cell0, cell_population);

        c_vector<double, 3> parent_location_tip_cell;
        parent_location_tip_cell = positions_tip_cell.first;
        TS_ASSERT_DELTA(norm_2(parent_location_tip_cell), norm_2(expected_parent_location_tip_cell), 1e-1);

        c_vector<double, 3> daughter_location_tip_cell;
        daughter_location_tip_cell = positions_tip_cell.second;
        TS_ASSERT_DELTA(norm_2(daughter_location_tip_cell), norm_2(expected_daughter_location_tip_cell), 1e-1);

        // archiving
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "sprouting_rule.arch";
        {
            SproutingRule* const p_sprouting_rule = new SproutingRule(0.08, 2.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_sprouting_rule;

            delete p_sprouting_rule;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            SproutingRule* p_sprouting_rule;
            input_arch >> p_sprouting_rule;

            TS_ASSERT_EQUALS(dynamic_cast<SproutingRule*>(p_sprouting_rule)->GetMaxSproutingRate(), 0.08);

            delete p_sprouting_rule;
        }
    }

    void NoTestSproutingRuleWithConstantVegfIn1D()
    {
        // creation of the mesh
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0u, false, 10, 5));
        nodes.push_back(new Node<1>(1u, false, 9, 5));
        nodes.push_back(new Node<1>(2u, false, 8, 5));

        NodesOnlyMesh<1> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<1> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        cell_population.Update();

        // Set up cell data on the cell population : initialisation for the division 
        for (typename NodeBasedCellPopulation<1>::Iterator cell_it = cell_population.Begin(); 
                                                cell_it != cell_population.End(); ++cell_it )
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_it);
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            c_vector<double, 1> node_coordinates = cell_population.GetLocationOfCellCentre(p_cell);

            double tortuosity = norm_2(node_coordinates);

            double t = SimulationTime::Instance()->GetTime(); 

            (*cell_it)->GetCellData()->SetItem("node_index", node_index);
            (*cell_it)->GetCellData()->SetItem("daughter_type", 0.0); // 0.0 = no daughter cell, 1.0 = daughter cell
            (*cell_it)->GetCellData()->SetItem("BranchNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingPoint", 0.0); // 0.0 = no, 1.0 = yes : necessary to modify the mutation state after sprouting in this modifier 
            (*cell_it)->GetCellData()->SetItem("DivisionNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingLeader", 0.0); // indice of the branching cell leading the branch of the cell considered 
            (*cell_it)->GetCellData()->SetItem("BranchingSegment", 0.0);
            (*cell_it)->GetCellData()->SetItem("ConsecutiveBranches", 1.0); 
            (*cell_it)->GetCellData()->SetItem("BirthTimeCell", t);
            (*cell_it)->GetCellData()->SetItem("TortuosityParent", tortuosity);
            (*cell_it)->GetCellData()->SetItem("TortuosityDaughter", tortuosity);
            (*cell_it)->GetCellData()->SetItem("LoopNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("Anastomosis", 0);
        }

        // create a vector of the expected coordinates of the parent and daughter cells : in this case, the tip cell divides into a vessel segment
        CellPtr p_cell0 = cell_population.GetCellUsingLocationIndex(2);
        c_vector<double,1> old_parent_location = cell_population.GetLocationOfCellCentre(p_cell0);

        // 0.3 is the MeinekeDivisionSeparation value 
        c_vector<double, 1> expected_daughter_location_vessel_segment  = old_parent_location;
        c_vector<double, 1> expected_parent_location_vessel_segment;
        expected_parent_location_vessel_segment[0] = old_parent_location[0] - 0.3*1;  

        c_vector<double, 1> expected_daughter_location_tip_cell;
        expected_daughter_location_tip_cell[0] = old_parent_location[0] - 0.5*0.3*(-1) + 0.5*0.3*1; 
        c_vector<double, 1> expected_parent_location_tip_cell;
        expected_parent_location_tip_cell[0] = old_parent_location[0] + 0.5*0.3*(-1) + 0.5*0.3*1;  

        // we create the sprouting rule object with a MaxSproutingRateConstant of 1.0
        typedef SproutingRuleWithConstantVegf<1,1> SproutingRuleWithConstantVegf;

        // for sourceterm = 0.1
        MAKE_PTR_ARGS(SproutingRuleWithConstantVegf, p_sprouting_rule_vessel_segment, (1.0, 2.0, 0.1, 1, 0.3, 0.98, 0.4));
        TS_ASSERT_DELTA(p_sprouting_rule_vessel_segment->GetSproutingProbabilityWithConstantVegf(), 0.013, 1e-3);

        // for sourceterm = 1.0
        MAKE_PTR_ARGS(SproutingRuleWithConstantVegf, p_sprouting_rule_tip_cell, (1.0, 2.0, 1.0, 1, 0.3, 0.98, 0.4));
        TS_ASSERT_DELTA(p_sprouting_rule_tip_cell->GetSproutingProbabilityWithConstantVegf(), 0.98, 1e-2);

        // case of vessel segment 
        // set the division rule to be the division rule for the cell population
        cell_population.SetCentreBasedDivisionRule(p_sprouting_rule_vessel_segment);
        boost::shared_ptr<AbstractCentreBasedDivisionRule<1,1> > p_sprouting_rule_vessel_segment_test = cell_population.GetCentreBasedDivisionRule();

        // test to check that the division rule returns the correct pair of vectors 
        std::pair<c_vector<double, 1>, c_vector<double, 1> > positions_vessel_segment = p_sprouting_rule_vessel_segment_test->CalculateCellDivisionVector(p_cell0, cell_population);

        c_vector<double, 1> parent_location_vessel_segment;
        parent_location_vessel_segment = positions_vessel_segment.first;
        TS_ASSERT_DELTA(parent_location_vessel_segment[0], expected_parent_location_vessel_segment[0], 1e-6);

        c_vector<double, 1> daughter_location_vessel_segment;
        daughter_location_vessel_segment = positions_vessel_segment.second;
        TS_ASSERT_DELTA(daughter_location_vessel_segment[0], expected_daughter_location_vessel_segment[0], 1e-6);

        // case of tip cell (take into account new position of the tip cell) 
        // set the division rule to be the division rule for the cell population
        cell_population.SetCentreBasedDivisionRule(p_sprouting_rule_tip_cell);
        boost::shared_ptr<AbstractCentreBasedDivisionRule<1,1> > p_sprouting_rule_tip_cell_test = cell_population.GetCentreBasedDivisionRule();

        // test to check that the division rule returns the correct pair of vectors 
        std::pair<c_vector<double, 1>, c_vector<double, 1> > positions_tip_cell = p_sprouting_rule_tip_cell_test->CalculateCellDivisionVector(p_cell0, cell_population);

        c_vector<double, 1> parent_location_tip_cell;
        parent_location_tip_cell = positions_tip_cell.first;
        TS_ASSERT_DELTA(parent_location_tip_cell[0], expected_parent_location_tip_cell[0], 1e-6);

        c_vector<double, 1> daughter_location_tip_cell;
        daughter_location_tip_cell = positions_tip_cell.second;
        TS_ASSERT_DELTA(daughter_location_tip_cell[0], expected_daughter_location_tip_cell[0], 1e-6);

        // archiving
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "sprouting_rule.arch";
        {
            SproutingRuleWithConstantVegf* const p_sprouting_rule = new SproutingRuleWithConstantVegf(0.08, 2.0, 0.5, 1, 0.3, 0.98, 0.4);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_sprouting_rule;

            delete p_sprouting_rule;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            SproutingRuleWithConstantVegf* p_sprouting_rule;
            input_arch >> p_sprouting_rule;

            TS_ASSERT_DELTA(dynamic_cast<SproutingRuleWithConstantVegf*>(p_sprouting_rule)->GetSproutingProbabilityWithConstantVegf(), 0.064, 1e-3);

            delete p_sprouting_rule;
        }
    }

    void NoTestSproutingRuleWithConstantVegfIn2D()
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
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

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

        cell_population.Update();

        // Set up cell data on the cell population : initialisation for the division 
        for (typename NodeBasedCellPopulation<2>::Iterator cell_it = cell_population.Begin(); 
                                                cell_it != cell_population.End(); ++cell_it )
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_it);
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            c_vector<double, 2> node_coordinates = cell_population.GetLocationOfCellCentre(p_cell);

            double tortuosity = norm_2(node_coordinates);

            double t = SimulationTime::Instance()->GetTime(); 

            (*cell_it)->GetCellData()->SetItem("node_index", node_index);
            (*cell_it)->GetCellData()->SetItem("daughter_type", 0.0); // 0.0 = no daughter cell, 1.0 = daughter cell
            (*cell_it)->GetCellData()->SetItem("BranchNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingPoint", 0.0); // 0.0 = no, 1.0 = yes : necessary to modify the mutation state after sprouting in this modifier 
            (*cell_it)->GetCellData()->SetItem("DivisionNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingLeader", 0.0); // indice of the branching cell leading the branch of the cell considered 
            (*cell_it)->GetCellData()->SetItem("BranchingSegment", 0.0);
            (*cell_it)->GetCellData()->SetItem("ConsecutiveBranches", 1.0); 
            (*cell_it)->GetCellData()->SetItem("BirthTimeCell", t);
            (*cell_it)->GetCellData()->SetItem("TortuosityParent", tortuosity);
            (*cell_it)->GetCellData()->SetItem("TortuosityDaughter", tortuosity);
            (*cell_it)->GetCellData()->SetItem("LoopNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("Anastomosis", 0);
        }

        // create a vector of the expected coordinates of the parent and daughter cells : in this case, the tip cell divides into a vessel segment
        CellPtr p_cell0 = cell_population.GetCellUsingLocationIndex(2);
        c_vector<double,2> old_parent_location = cell_population.GetLocationOfCellCentre(p_cell0);

        // 0.3 is the MeinekeDivisionSeparation value and (-1, 0) is the daughter_direction
        c_vector<double, 2> expected_daughter_location_vessel_segment  = old_parent_location;
        c_vector<double, 2> expected_parent_location_vessel_segment;
        expected_parent_location_vessel_segment[0] = old_parent_location[0] - 0.3*1;  
        expected_parent_location_vessel_segment[1] = old_parent_location[1] - 0.3*0;

        c_vector<double, 2> expected_daughter_location_tip_cell;
        expected_daughter_location_tip_cell[0] = old_parent_location[0] - 0.5*0.3*0 + 0.5*0.3*1;  
        expected_daughter_location_tip_cell[1] = old_parent_location[1] - 0.5*0.3*1 + 0.5*0.3*0;
        c_vector<double, 2> expected_parent_location_tip_cell;
        expected_parent_location_tip_cell[0] = old_parent_location[0] + 0.5*0.3*0 + 0.5*0.3*1;  
        expected_parent_location_tip_cell[1] = old_parent_location[1] + 0.5*0.3*1 + 0.5*0.3*0;

        // we create the sprouting rule object with a MaxSproutingRateConstant of 1.0
        typedef SproutingRuleWithConstantVegf<2,2> SproutingRuleWithConstantVegf;

        // for sourceterm = 0.1
        MAKE_PTR_ARGS(SproutingRuleWithConstantVegf, p_sprouting_rule_vessel_segment, (1.0, 2.0, 0.1, 1, 0.3, 0.98, 0.4));
        TS_ASSERT_DELTA(p_sprouting_rule_vessel_segment->GetSproutingProbabilityWithConstantVegf(), 0.013, 1e-3);

        // for sourceterm = 1.0
        MAKE_PTR_ARGS(SproutingRuleWithConstantVegf, p_sprouting_rule_tip_cell, (1.0, 2.0, 1.0, 1, 0.3, 0.98, 0.4));
        TS_ASSERT_DELTA(p_sprouting_rule_tip_cell->GetSproutingProbabilityWithConstantVegf(), 0.98, 1e-2);

        // case of vessel segment 
        // set the division rule to be the division rule for the cell population
        cell_population.SetCentreBasedDivisionRule(p_sprouting_rule_vessel_segment);
        boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_sprouting_rule_vessel_segment_test = cell_population.GetCentreBasedDivisionRule();

        // test to check that the division rule returns the correct pair of vectors 
        std::pair<c_vector<double, 2>, c_vector<double, 2> > positions_vessel_segment = p_sprouting_rule_vessel_segment_test->CalculateCellDivisionVector(p_cell0, cell_population);

        c_vector<double, 2> parent_location_vessel_segment;
        parent_location_vessel_segment = positions_vessel_segment.first;
        TS_ASSERT_DELTA(parent_location_vessel_segment[0], expected_parent_location_vessel_segment[0], 1e-6);
        TS_ASSERT_DELTA(parent_location_vessel_segment[1], expected_parent_location_vessel_segment[1], 1e-6);

        c_vector<double, 2> daughter_location_vessel_segment;
        daughter_location_vessel_segment = positions_vessel_segment.second;
        TS_ASSERT_DELTA(daughter_location_vessel_segment[0], expected_daughter_location_vessel_segment[0], 1e-6);
        TS_ASSERT_DELTA(daughter_location_vessel_segment[1], expected_daughter_location_vessel_segment[1], 1e-6);

        // case of tip cell (take into account new position of the tip cell) 
        // set the division rule to be the division rule for the cell population
        cell_population.SetCentreBasedDivisionRule(p_sprouting_rule_tip_cell);
        boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_sprouting_rule_tip_cell_test = cell_population.GetCentreBasedDivisionRule();

        // test to check that the division rule returns the correct pair of vectors 
        std::pair<c_vector<double, 2>, c_vector<double, 2> > positions_tip_cell = p_sprouting_rule_tip_cell_test->CalculateCellDivisionVector(p_cell0, cell_population);

        c_vector<double, 2> parent_location_tip_cell;
        parent_location_tip_cell = positions_tip_cell.first;
        TS_ASSERT_DELTA(parent_location_tip_cell[0], expected_parent_location_tip_cell[0], 1e-6);
        TS_ASSERT_DELTA(parent_location_tip_cell[1], expected_parent_location_tip_cell[1], 1e-6);

        c_vector<double, 2> daughter_location_tip_cell;
        daughter_location_tip_cell = positions_tip_cell.second;
        TS_ASSERT_DELTA(daughter_location_tip_cell[0], expected_daughter_location_tip_cell[0], 1e-6);
        TS_ASSERT_DELTA(daughter_location_tip_cell[1], expected_daughter_location_tip_cell[1], 1e-6);

        // archiving
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "sprouting_rule.arch";
        {
            SproutingRuleWithConstantVegf* const p_sprouting_rule = new SproutingRuleWithConstantVegf(0.08, 2.0, 0.5, 1, 0.3, 0.98, 0.4);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_sprouting_rule;

            delete p_sprouting_rule;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            SproutingRuleWithConstantVegf* p_sprouting_rule;
            input_arch >> p_sprouting_rule;

            TS_ASSERT_DELTA(dynamic_cast<SproutingRuleWithConstantVegf*>(p_sprouting_rule)->GetSproutingProbabilityWithConstantVegf(), 0.064, 1e-3);

            delete p_sprouting_rule;
        }
    }

    void NoTestSproutingRuleWithConstantVegfIn3D()
    {
        // creation of the mesh
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u, false, 10, 5));
        nodes.push_back(new Node<3>(1u, false, 9, 5));
        nodes.push_back(new Node<3>(2u, false, 8, 5));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        cell_population.Update();

        // Set up cell data on the cell population : initialisation for the division 
        for (typename NodeBasedCellPopulation<3>::Iterator cell_it = cell_population.Begin(); 
                                                cell_it != cell_population.End(); ++cell_it )
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_it);
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            c_vector<double, 3> node_coordinates = cell_population.GetLocationOfCellCentre(p_cell);

            double tortuosity = norm_2(node_coordinates);

            double t = SimulationTime::Instance()->GetTime(); 

            (*cell_it)->GetCellData()->SetItem("node_index", node_index);
            (*cell_it)->GetCellData()->SetItem("daughter_type", 0.0); // 0.0 = no daughter cell, 1.0 = daughter cell
            (*cell_it)->GetCellData()->SetItem("BranchNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingPoint", 0.0); // 0.0 = no, 1.0 = yes : necessary to modify the mutation state after sprouting in this modifier 
            (*cell_it)->GetCellData()->SetItem("DivisionNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingLeader", 0.0); // indice of the branching cell leading the branch of the cell considered 
            (*cell_it)->GetCellData()->SetItem("BranchingSegment", 0.0);
            (*cell_it)->GetCellData()->SetItem("ConsecutiveBranches", 1.0); 
            (*cell_it)->GetCellData()->SetItem("BirthTimeCell", t);
            (*cell_it)->GetCellData()->SetItem("TortuosityParent", tortuosity);
            (*cell_it)->GetCellData()->SetItem("TortuosityDaughter", tortuosity);
            (*cell_it)->GetCellData()->SetItem("LoopNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("Anastomosis", 0);
        }

        // create a vector of the expected coordinates of the parent and daughter cells : in this case, the tip cell divides into a vessel segment
        CellPtr p_cell0 = cell_population.GetCellUsingLocationIndex(2);
        c_vector<double,3> old_parent_location = cell_population.GetLocationOfCellCentre(p_cell0);

        // 0.3 is the MeinekeDivisionSeparation value and (-1, 0) is the daughter_direction
        c_vector<double, 3> expected_daughter_location_vessel_segment  = old_parent_location;
        c_vector<double, 3> expected_parent_location_vessel_segment;
        expected_parent_location_vessel_segment[0] = old_parent_location[0] - 0.3*1;  
        expected_parent_location_vessel_segment[1] = old_parent_location[1] - 0.3*0;
        expected_parent_location_vessel_segment[2] = old_parent_location[2] - 0.3*0;

        c_vector<double, 3> expected_daughter_location_tip_cell;
        expected_daughter_location_tip_cell[0] = old_parent_location[0] - 0.5*0.3*0 + 0.5*0.3*1;  
        expected_daughter_location_tip_cell[1] = old_parent_location[1] - 0.5*0.3*1 + 0.5*0.3*0;
        expected_daughter_location_tip_cell[2] = old_parent_location[2] - 0.5*0.3*0 + 0.5*0.3*0;
        c_vector<double, 3> expected_parent_location_tip_cell;
        expected_parent_location_tip_cell[0] = old_parent_location[0] + 0.5*0.3*0 + 0.5*0.3*1;  
        expected_parent_location_tip_cell[1] = old_parent_location[1] + 0.5*0.3*1 + 0.5*0.3*0;
        expected_parent_location_tip_cell[2] = old_parent_location[2] + 0.5*0.3*0 + 0.5*0.3*0;

        // we create the sprouting rule object with a MaxSproutingRateConstant of 1.0
        typedef SproutingRuleWithConstantVegf<3,3> SproutingRuleWithConstantVegf;

        // for sourceterm = 0.1
        MAKE_PTR_ARGS(SproutingRuleWithConstantVegf, p_sprouting_rule_vessel_segment, (1.0, 2.0, 0.1, 1, 0.3, 0.98, 0.4));
        TS_ASSERT_DELTA(p_sprouting_rule_vessel_segment->GetSproutingProbabilityWithConstantVegf(), 0.013, 1e-3);

        // for sourceterm = 1.0
        MAKE_PTR_ARGS(SproutingRuleWithConstantVegf, p_sprouting_rule_tip_cell, (1.0, 2.0, 1.0, 1, 0.3, 0.98, 0.4));
        TS_ASSERT_DELTA(p_sprouting_rule_tip_cell->GetSproutingProbabilityWithConstantVegf(), 0.98, 1e-2);

        // case of vessel segment 
        // set the division rule to be the division rule for the cell population
        cell_population.SetCentreBasedDivisionRule(p_sprouting_rule_vessel_segment);
        boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_sprouting_rule_vessel_segment_test = cell_population.GetCentreBasedDivisionRule();

        // test to check that the division rule returns the correct pair of vectors 
        std::pair<c_vector<double, 3>, c_vector<double, 3> > positions_vessel_segment = p_sprouting_rule_vessel_segment_test->CalculateCellDivisionVector(p_cell0, cell_population);

        c_vector<double, 3> parent_location_vessel_segment;
        parent_location_vessel_segment = positions_vessel_segment.first;
        TS_ASSERT_DELTA(parent_location_vessel_segment[0], expected_parent_location_vessel_segment[0], 1e-6);
        TS_ASSERT_DELTA(parent_location_vessel_segment[1], expected_parent_location_vessel_segment[1], 1e-6);
        TS_ASSERT_DELTA(parent_location_vessel_segment[2], expected_parent_location_vessel_segment[2], 1e-6);

        c_vector<double, 3> daughter_location_vessel_segment;
        daughter_location_vessel_segment = positions_vessel_segment.second;
        TS_ASSERT_DELTA(daughter_location_vessel_segment[0], expected_daughter_location_vessel_segment[0], 1e-6);
        TS_ASSERT_DELTA(daughter_location_vessel_segment[1], expected_daughter_location_vessel_segment[1], 1e-6);
        TS_ASSERT_DELTA(daughter_location_vessel_segment[2], expected_daughter_location_vessel_segment[2], 1e-6);

        // case of tip cell (take into account new position of the tip cell) 
        // set the division rule to be the division rule for the cell population
        cell_population.SetCentreBasedDivisionRule(p_sprouting_rule_tip_cell);
        boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_sprouting_rule_tip_cell_test = cell_population.GetCentreBasedDivisionRule();

        // test to check that the division rule returns the correct pair of vectors 
        std::pair<c_vector<double, 3>, c_vector<double, 3> > positions_tip_cell = p_sprouting_rule_tip_cell_test->CalculateCellDivisionVector(p_cell0, cell_population);

        c_vector<double, 3> parent_location_tip_cell;
        parent_location_tip_cell = positions_tip_cell.first;
        TS_ASSERT_DELTA(norm_2(parent_location_tip_cell), norm_2(expected_parent_location_tip_cell), 1e-1);

        c_vector<double, 3> daughter_location_tip_cell;
        daughter_location_tip_cell = positions_tip_cell.second;
        TS_ASSERT_DELTA(norm_2(daughter_location_tip_cell), norm_2(expected_daughter_location_tip_cell), 1e-1);

        // archiving
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "sprouting_rule.arch";
        {
            SproutingRuleWithConstantVegf* const p_sprouting_rule = new SproutingRuleWithConstantVegf(0.08, 2.0, 0.5, 1, 0.3, 0.98, 0.4);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_sprouting_rule;

            delete p_sprouting_rule;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            SproutingRuleWithConstantVegf* p_sprouting_rule;
            input_arch >> p_sprouting_rule;

            TS_ASSERT_DELTA(dynamic_cast<SproutingRuleWithConstantVegf*>(p_sprouting_rule)->GetSproutingProbabilityWithConstantVegf(), 0.064, 1e-3);

            delete p_sprouting_rule;
        }
    }

    void NoTestSproutingRuleWithAnalyticalApproximationPdeIn1D()
    {
        // creation of the mesh
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0u, false, 10, 5));
        nodes.push_back(new Node<1>(1u, false, 9, 5));
        nodes.push_back(new Node<1>(2u, false, 8, 5));

        NodesOnlyMesh<1> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<1> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        cell_population.Update();

        // Set up cell data on the cell population : initialisation for the division 
        for (typename NodeBasedCellPopulation<1>::Iterator cell_it = cell_population.Begin(); 
                                                cell_it != cell_population.End(); ++cell_it )
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_it);
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            c_vector<double, 1> node_coordinates = cell_population.GetLocationOfCellCentre(p_cell);

            double tortuosity = norm_2(node_coordinates);

            double t = SimulationTime::Instance()->GetTime(); 

            (*cell_it)->GetCellData()->SetItem("node_index", node_index);
            (*cell_it)->GetCellData()->SetItem("daughter_type", 0.0); // 0.0 = no daughter cell, 1.0 = daughter cell
            (*cell_it)->GetCellData()->SetItem("BranchNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingPoint", 0.0); // 0.0 = no, 1.0 = yes : necessary to modify the mutation state after sprouting in this modifier 
            (*cell_it)->GetCellData()->SetItem("DivisionNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingLeader", 0.0); // indice of the branching cell leading the branch of the cell considered 
            (*cell_it)->GetCellData()->SetItem("BranchingSegment", 0.0);
            (*cell_it)->GetCellData()->SetItem("ConsecutiveBranches", 1.0); 
            (*cell_it)->GetCellData()->SetItem("BirthTimeCell", t);
            (*cell_it)->GetCellData()->SetItem("TortuosityParent", tortuosity);
            (*cell_it)->GetCellData()->SetItem("TortuosityDaughter", tortuosity);
            (*cell_it)->GetCellData()->SetItem("LoopNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("Anastomosis", 0);
        }

        // create a vector of the expected coordinates of the parent and daughter cells : in this case, the tip cell divides into a vessel segment
        CellPtr p_cell0 = cell_population.GetCellUsingLocationIndex(2);
        c_vector<double,1> old_parent_location = cell_population.GetLocationOfCellCentre(p_cell0);

        // 0.3 is the MeinekeDivisionSeparation value and (-1, 0) is the daughter_direction
        c_vector<double, 1> expected_daughter_location_vessel_segment  = old_parent_location;
        c_vector<double, 1> expected_parent_location_vessel_segment;
        expected_parent_location_vessel_segment[0] = old_parent_location[0] - 0.3*1;  

        c_vector<double, 1> expected_daughter_location_tip_cell;
        expected_daughter_location_tip_cell[0] = old_parent_location[0] - 0.5*0.3*(-1) + 0.5*0.3*1;  
        c_vector<double, 1> expected_parent_location_tip_cell;
        expected_parent_location_tip_cell[0] = old_parent_location[0] + 0.5*0.3*(-1) + 0.5*0.3*1;  

        // we create the sprouting rule object with a MaxSproutingRateConstant of 1.0
        typedef SproutingRuleWithAnalyticalApproximationPde<1,1> SproutingRuleWithAnalyticalApproximationPde;

        // for sourceterm = 0.1
        MAKE_PTR_ARGS(SproutingRuleWithAnalyticalApproximationPde, p_sprouting_rule_vessel_segment, (1.0, 2.0, 1e4, 1.0, 0.0, 1.0, 0.1, 0.1, 1, 0.3, 0.98, 0.4));
        TS_ASSERT_DELTA(p_sprouting_rule_vessel_segment->GetSproutingProbability(cell_population, p_cell0), 0.013, 1e-4);

        // for sourceterm = 1.0
        MAKE_PTR_ARGS(SproutingRuleWithAnalyticalApproximationPde, p_sprouting_rule_tip_cell, (1.0, 2.0, 1e4, 1.0, 0.0, 1.0, 1.0, 0.1, 1, 0.3, 0.98, 0.4));
        TS_ASSERT_DELTA(p_sprouting_rule_tip_cell->GetSproutingProbability(cell_population, p_cell0), 0.97, 1e-2);

        // case of vessel segment 
        // set the division rule to be the division rule for the cell population
        cell_population.SetCentreBasedDivisionRule(p_sprouting_rule_vessel_segment);
        boost::shared_ptr<AbstractCentreBasedDivisionRule<1,1> > p_sprouting_rule_vessel_segment_test = cell_population.GetCentreBasedDivisionRule();

        // test to check that the division rule returns the correct pair of vectors 
        std::pair<c_vector<double, 1>, c_vector<double, 1> > positions_vessel_segment = p_sprouting_rule_vessel_segment_test->CalculateCellDivisionVector(p_cell0, cell_population);

        c_vector<double, 1> parent_location_vessel_segment;
        parent_location_vessel_segment = positions_vessel_segment.first;
        TS_ASSERT_DELTA(parent_location_vessel_segment[0], expected_parent_location_vessel_segment[0], 1e-6);

        c_vector<double, 1> daughter_location_vessel_segment;
        daughter_location_vessel_segment = positions_vessel_segment.second;
        TS_ASSERT_DELTA(daughter_location_vessel_segment[0], expected_daughter_location_vessel_segment[0], 1e-6);

        // case of tip cell (take into account new position of the tip cell) 
        // set the division rule to be the division rule for the cell population
        cell_population.SetCentreBasedDivisionRule(p_sprouting_rule_tip_cell);
        boost::shared_ptr<AbstractCentreBasedDivisionRule<1,1> > p_sprouting_rule_tip_cell_test = cell_population.GetCentreBasedDivisionRule();

        // test to check that the division rule returns the correct pair of vectors 
        std::pair<c_vector<double, 1>, c_vector<double, 1> > positions_tip_cell = p_sprouting_rule_tip_cell_test->CalculateCellDivisionVector(p_cell0, cell_population);

        c_vector<double, 1> parent_location_tip_cell;
        parent_location_tip_cell = positions_tip_cell.first;
        TS_ASSERT_DELTA(parent_location_tip_cell[0], expected_parent_location_tip_cell[0], 1e-6);

        c_vector<double, 1> daughter_location_tip_cell;
        daughter_location_tip_cell = positions_tip_cell.second;
        TS_ASSERT_DELTA(daughter_location_tip_cell[0], expected_daughter_location_tip_cell[0], 1e-6);

        // archiving
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "sprouting_rule.arch";
        {
            SproutingRuleWithAnalyticalApproximationPde* const p_sprouting_rule = new SproutingRuleWithAnalyticalApproximationPde(0.08, 2.0, 1e4, 1.0, 0.0, 1.0, 0.5, 0.1, 1, 0.3, 0.98, 0.4);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_sprouting_rule;

            delete p_sprouting_rule;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            SproutingRuleWithAnalyticalApproximationPde* p_sprouting_rule;
            input_arch >> p_sprouting_rule;

            delete p_sprouting_rule;
        }
    }

    void NoTestSproutingRuleWithAnalyticalApproximationPdeIn2D()
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
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

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

        cell_population.Update();

        // Set up cell data on the cell population : initialisation for the division 
        for (typename NodeBasedCellPopulation<2>::Iterator cell_it = cell_population.Begin(); 
                                                cell_it != cell_population.End(); ++cell_it )
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_it);
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            c_vector<double, 2> node_coordinates = cell_population.GetLocationOfCellCentre(p_cell);

            double tortuosity = norm_2(node_coordinates);

            double t = SimulationTime::Instance()->GetTime(); 

            (*cell_it)->GetCellData()->SetItem("node_index", node_index);
            (*cell_it)->GetCellData()->SetItem("daughter_type", 0.0); // 0.0 = no daughter cell, 1.0 = daughter cell
            (*cell_it)->GetCellData()->SetItem("BranchNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingPoint", 0.0); // 0.0 = no, 1.0 = yes : necessary to modify the mutation state after sprouting in this modifier 
            (*cell_it)->GetCellData()->SetItem("DivisionNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingLeader", 0.0); // indice of the branching cell leading the branch of the cell considered 
            (*cell_it)->GetCellData()->SetItem("BranchingSegment", 0.0);
            (*cell_it)->GetCellData()->SetItem("ConsecutiveBranches", 1.0); 
            (*cell_it)->GetCellData()->SetItem("BirthTimeCell", t);
            (*cell_it)->GetCellData()->SetItem("TortuosityParent", tortuosity);
            (*cell_it)->GetCellData()->SetItem("TortuosityDaughter", tortuosity);
            (*cell_it)->GetCellData()->SetItem("LoopNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("Anastomosis", 0);
        }

        // create a vector of the expected coordinates of the parent and daughter cells : in this case, the tip cell divides into a vessel segment
        CellPtr p_cell0 = cell_population.GetCellUsingLocationIndex(2);
        c_vector<double,2> old_parent_location = cell_population.GetLocationOfCellCentre(p_cell0);

        // 0.3 is the MeinekeDivisionSeparation value and (-1, 0) is the daughter_direction
        c_vector<double, 2> expected_daughter_location_vessel_segment  = old_parent_location;
        c_vector<double, 2> expected_parent_location_vessel_segment;
        expected_parent_location_vessel_segment[0] = old_parent_location[0] - 0.3*1;  
        expected_parent_location_vessel_segment[1] = old_parent_location[1] - 0.3*0;

        c_vector<double, 2> expected_daughter_location_tip_cell;
        expected_daughter_location_tip_cell[0] = old_parent_location[0] - 0.5*0.3*0 + 0.5*0.3*1;  
        expected_daughter_location_tip_cell[1] = old_parent_location[1] - 0.5*0.3*1 + 0.5*0.3*0;
        c_vector<double, 2> expected_parent_location_tip_cell;
        expected_parent_location_tip_cell[0] = old_parent_location[0] + 0.5*0.3*0 + 0.5*0.3*1;  
        expected_parent_location_tip_cell[1] = old_parent_location[1] + 0.5*0.3*1 + 0.5*0.3*0;

        // we create the sprouting rule object with a MaxSproutingRateConstant of 1.0
        typedef SproutingRuleWithAnalyticalApproximationPde<2,2> SproutingRuleWithAnalyticalApproximationPde;

        // for sourceterm = 0.1
        MAKE_PTR_ARGS(SproutingRuleWithAnalyticalApproximationPde, p_sprouting_rule_vessel_segment, (1.0, 2.0, 1e4, 1.0, 0.0, 1.0, 0.1, 0.1, 1, 0.3, 0.98, 0.4));
        TS_ASSERT_DELTA(p_sprouting_rule_vessel_segment->GetSproutingProbability(cell_population, p_cell0), 0.013, 1e-4);

        // for sourceterm = 1.0
        MAKE_PTR_ARGS(SproutingRuleWithAnalyticalApproximationPde, p_sprouting_rule_tip_cell, (1.0, 2.0, 1e4, 1.0, 0.0, 1.0, 1.0, 0.1, 1, 0.3, 0.98, 0.4));
        TS_ASSERT_DELTA(p_sprouting_rule_tip_cell->GetSproutingProbability(cell_population, p_cell0), 0.97, 1e-2);

        // case of vessel segment 
        // set the division rule to be the division rule for the cell population
        cell_population.SetCentreBasedDivisionRule(p_sprouting_rule_vessel_segment);
        boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_sprouting_rule_vessel_segment_test = cell_population.GetCentreBasedDivisionRule();

        // test to check that the division rule returns the correct pair of vectors 
        std::pair<c_vector<double, 2>, c_vector<double, 2> > positions_vessel_segment = p_sprouting_rule_vessel_segment_test->CalculateCellDivisionVector(p_cell0, cell_population);

        c_vector<double, 2> parent_location_vessel_segment;
        parent_location_vessel_segment = positions_vessel_segment.first;
        TS_ASSERT_DELTA(parent_location_vessel_segment[0], expected_parent_location_vessel_segment[0], 1e-6);
        TS_ASSERT_DELTA(parent_location_vessel_segment[1], expected_parent_location_vessel_segment[1], 1e-6);

        c_vector<double, 2> daughter_location_vessel_segment;
        daughter_location_vessel_segment = positions_vessel_segment.second;
        TS_ASSERT_DELTA(daughter_location_vessel_segment[0], expected_daughter_location_vessel_segment[0], 1e-6);
        TS_ASSERT_DELTA(daughter_location_vessel_segment[1], expected_daughter_location_vessel_segment[1], 1e-6);

        // case of tip cell (take into account new position of the tip cell) 
        // set the division rule to be the division rule for the cell population
        cell_population.SetCentreBasedDivisionRule(p_sprouting_rule_tip_cell);
        boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_sprouting_rule_tip_cell_test = cell_population.GetCentreBasedDivisionRule();

        // test to check that the division rule returns the correct pair of vectors 
        std::pair<c_vector<double, 2>, c_vector<double, 2> > positions_tip_cell = p_sprouting_rule_tip_cell_test->CalculateCellDivisionVector(p_cell0, cell_population);

        c_vector<double, 2> parent_location_tip_cell;
        parent_location_tip_cell = positions_tip_cell.first;
        TS_ASSERT_DELTA(parent_location_tip_cell[0], expected_parent_location_tip_cell[0], 1e-6);
        TS_ASSERT_DELTA(parent_location_tip_cell[1], expected_parent_location_tip_cell[1], 1e-6);

        c_vector<double, 2> daughter_location_tip_cell;
        daughter_location_tip_cell = positions_tip_cell.second;
        TS_ASSERT_DELTA(daughter_location_tip_cell[0], expected_daughter_location_tip_cell[0], 1e-6);
        TS_ASSERT_DELTA(daughter_location_tip_cell[1], expected_daughter_location_tip_cell[1], 1e-6);

        // archiving
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "sprouting_rule.arch";
        {
            SproutingRuleWithAnalyticalApproximationPde* const p_sprouting_rule = new SproutingRuleWithAnalyticalApproximationPde(0.08, 2.0, 1e4, 1.0, 0.0, 1.0, 0.5, 0.1, 1, 0.3, 0.98, 0.4);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_sprouting_rule;

            delete p_sprouting_rule;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            SproutingRuleWithAnalyticalApproximationPde* p_sprouting_rule;
            input_arch >> p_sprouting_rule;

            delete p_sprouting_rule;
        }
    }

    void NoTestSproutingRuleWithAnalyticalApproximationPdeIn3D()
    {
        // creation of the mesh
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u, false, 10, 5));
        nodes.push_back(new Node<3>(1u, false, 9, 5));
        nodes.push_back(new Node<3>(2u, false, 8, 5));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        cell_population.Update();

        // Set up cell data on the cell population : initialisation for the division 
        for (typename NodeBasedCellPopulation<3>::Iterator cell_it = cell_population.Begin(); 
                                                cell_it != cell_population.End(); ++cell_it )
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_it);
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            c_vector<double, 3> node_coordinates = cell_population.GetLocationOfCellCentre(p_cell);

            double tortuosity = norm_2(node_coordinates);

            double t = SimulationTime::Instance()->GetTime(); 

            (*cell_it)->GetCellData()->SetItem("node_index", node_index);
            (*cell_it)->GetCellData()->SetItem("daughter_type", 0.0); // 0.0 = no daughter cell, 1.0 = daughter cell
            (*cell_it)->GetCellData()->SetItem("BranchNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingPoint", 0.0); // 0.0 = no, 1.0 = yes : necessary to modify the mutation state after sprouting in this modifier 
            (*cell_it)->GetCellData()->SetItem("DivisionNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingLeader", 0.0); // indice of the branching cell leading the branch of the cell considered 
            (*cell_it)->GetCellData()->SetItem("BranchingSegment", 0.0);
            (*cell_it)->GetCellData()->SetItem("ConsecutiveBranches", 1.0); 
            (*cell_it)->GetCellData()->SetItem("BirthTimeCell", t);
            (*cell_it)->GetCellData()->SetItem("TortuosityParent", tortuosity);
            (*cell_it)->GetCellData()->SetItem("TortuosityDaughter", tortuosity);
            (*cell_it)->GetCellData()->SetItem("LoopNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("Anastomosis", 0);
        }

        // create a vector of the expected coordinates of the parent and daughter cells : in this case, the tip cell divides into a vessel segment
        CellPtr p_cell0 = cell_population.GetCellUsingLocationIndex(2);
        c_vector<double,3> old_parent_location = cell_population.GetLocationOfCellCentre(p_cell0);

        // 0.3 is the MeinekeDivisionSeparation value 
        c_vector<double, 3> expected_daughter_location_vessel_segment  = old_parent_location;
        c_vector<double, 3> expected_parent_location_vessel_segment;
        expected_parent_location_vessel_segment[0] = old_parent_location[0] - 0.3*1;  
        expected_parent_location_vessel_segment[1] = old_parent_location[1] - 0.3*0;
        expected_parent_location_vessel_segment[2] = old_parent_location[2] - 0.3*0;

        c_vector<double, 3> expected_daughter_location_tip_cell;
        expected_daughter_location_tip_cell[0] = old_parent_location[0] - 0.5*0.3*0 + 0.5*0.3*1;  
        expected_daughter_location_tip_cell[1] = old_parent_location[1] - 0.5*0.3*1 + 0.5*0.3*0;
        expected_daughter_location_tip_cell[2] = old_parent_location[2] - 0.5*0.3*0 + 0.5*0.3*0;
        c_vector<double, 3> expected_parent_location_tip_cell;
        expected_parent_location_tip_cell[0] = old_parent_location[0] + 0.5*0.3*0 + 0.5*0.3*1;  
        expected_parent_location_tip_cell[1] = old_parent_location[1] + 0.5*0.3*1 + 0.5*0.3*0;
        expected_parent_location_tip_cell[2] = old_parent_location[2] + 0.5*0.3*0 + 0.5*0.3*0;

        // we create the sprouting rule object with a MaxSproutingRateConstant of 1.0
        typedef SproutingRuleWithAnalyticalApproximationPde<3,3> SproutingRuleWithAnalyticalApproximationPde;

        // for sourceterm = 0.1
        MAKE_PTR_ARGS(SproutingRuleWithAnalyticalApproximationPde, p_sprouting_rule_vessel_segment, (1.0, 2.0, 1e4, 1.0, 0.0, 1.0, 0.1, 0.1, 1, 0.3, 0.98, 0.4));
        TS_ASSERT_DELTA(p_sprouting_rule_vessel_segment->GetSproutingProbability(cell_population, p_cell0), 0.013, 1e-4);

        // for sourceterm = 1.0
        MAKE_PTR_ARGS(SproutingRuleWithAnalyticalApproximationPde, p_sprouting_rule_tip_cell, (1.0, 2.0, 1e4, 1.0, 0.0, 1.0, 1.0, 0.1, 1, 0.3, 0.98, 0.4));
        TS_ASSERT_DELTA(p_sprouting_rule_tip_cell->GetSproutingProbability(cell_population, p_cell0), 0.97, 1e-2);

        // case of vessel segment 
        // set the division rule to be the division rule for the cell population
        cell_population.SetCentreBasedDivisionRule(p_sprouting_rule_vessel_segment);
        boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_sprouting_rule_vessel_segment_test = cell_population.GetCentreBasedDivisionRule();

        // test to check that the division rule returns the correct pair of vectors 
        std::pair<c_vector<double, 3>, c_vector<double, 3> > positions_vessel_segment = p_sprouting_rule_vessel_segment_test->CalculateCellDivisionVector(p_cell0, cell_population);

        c_vector<double, 3> parent_location_vessel_segment;
        parent_location_vessel_segment = positions_vessel_segment.first;
        TS_ASSERT_DELTA(parent_location_vessel_segment[0], expected_parent_location_vessel_segment[0], 1e-6);
        TS_ASSERT_DELTA(parent_location_vessel_segment[1], expected_parent_location_vessel_segment[1], 1e-6);
        TS_ASSERT_DELTA(parent_location_vessel_segment[2], expected_parent_location_vessel_segment[2], 1e-6);

        c_vector<double, 3> daughter_location_vessel_segment;
        daughter_location_vessel_segment = positions_vessel_segment.second;
        TS_ASSERT_DELTA(daughter_location_vessel_segment[0], expected_daughter_location_vessel_segment[0], 1e-6);
        TS_ASSERT_DELTA(daughter_location_vessel_segment[1], expected_daughter_location_vessel_segment[1], 1e-6);
        TS_ASSERT_DELTA(daughter_location_vessel_segment[2], expected_daughter_location_vessel_segment[2], 1e-6);

        // case of tip cell (take into account new position of the tip cell) 
        // set the division rule to be the division rule for the cell population
        cell_population.SetCentreBasedDivisionRule(p_sprouting_rule_tip_cell);
        boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_sprouting_rule_tip_cell_test = cell_population.GetCentreBasedDivisionRule();

        // test to check that the division rule returns the correct pair of vectors 
        std::pair<c_vector<double, 3>, c_vector<double, 3> > positions_tip_cell = p_sprouting_rule_tip_cell_test->CalculateCellDivisionVector(p_cell0, cell_population);

        c_vector<double, 3> parent_location_tip_cell;
        parent_location_tip_cell = positions_tip_cell.first;
        TS_ASSERT_DELTA(norm_2(parent_location_tip_cell), norm_2(expected_parent_location_tip_cell), 1e-1);

        c_vector<double, 3> daughter_location_tip_cell;
        daughter_location_tip_cell = positions_tip_cell.second;
        TS_ASSERT_DELTA(norm_2(daughter_location_tip_cell), norm_2(expected_daughter_location_tip_cell), 1e-1);

        // archiving
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "sprouting_rule.arch";
        {
            SproutingRuleWithAnalyticalApproximationPde* const p_sprouting_rule = new SproutingRuleWithAnalyticalApproximationPde(0.08, 2.0, 1e4, 1.0, 0.0, 1.0, 0.5, 0.1, 1, 0.3, 0.98, 0.4);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_sprouting_rule;

            delete p_sprouting_rule;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            SproutingRuleWithAnalyticalApproximationPde* p_sprouting_rule;
            input_arch >> p_sprouting_rule;

            delete p_sprouting_rule;
        }
    }

    void NoTestSproutingRule1WithPdesIn1D() 
    {
        // creation of the mesh
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0u, false, 10, 5, 5));
        nodes.push_back(new Node<1>(1u, false, 9, 5, 5));
        nodes.push_back(new Node<1>(2u, false, 8, 5, 5));

        NodesOnlyMesh<1> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<1> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        cell_population.Update();

        // Set up cell data on the cell population : initialisation for the division 
        for (typename NodeBasedCellPopulation<1>::Iterator cell_it = cell_population.Begin(); 
                                                cell_it != cell_population.End(); ++cell_it )
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_it);
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            c_vector<double, 1> node_coordinates = cell_population.GetLocationOfCellCentre(p_cell);

            double tortuosity = norm_2(node_coordinates);

            double t = SimulationTime::Instance()->GetTime(); 

            (*cell_it)->GetCellData()->SetItem("node_index", node_index);
            (*cell_it)->GetCellData()->SetItem("daughter_type", 0.0); // 0.0 = no daughter cell, 1.0 = daughter cell
            (*cell_it)->GetCellData()->SetItem("BranchNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingPoint", 0.0); // 0.0 = no, 1.0 = yes : necessary to modify the mutation state after sprouting in this modifier 
            (*cell_it)->GetCellData()->SetItem("DivisionNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingLeader", 0.0); // indice of the branching cell leading the branch of the cell considered 
            (*cell_it)->GetCellData()->SetItem("BranchingSegment", 0.0);
            (*cell_it)->GetCellData()->SetItem("ConsecutiveBranches", 1.0); 
            (*cell_it)->GetCellData()->SetItem("BirthTimeCell", t);
            (*cell_it)->GetCellData()->SetItem("TortuosityParent", tortuosity);
            (*cell_it)->GetCellData()->SetItem("TortuosityDaughter", tortuosity);
            (*cell_it)->GetCellData()->SetItem("LoopNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("Anastomosis", 0);
        }

        // create a vector of the expected coordinates of the parent and daughter cells : in this case, the tip cell divides into a vessel segment
        CellPtr p_cell_forvesselsegment = cell_population.GetCellUsingLocationIndex(2);
        c_vector<double,1> old_parent_location_forvesselsegment = cell_population.GetLocationOfCellCentre(p_cell_forvesselsegment);

        ChastePoint<1> lower(0.0);
        ChastePoint<1> upper(10.0);
        MAKE_PTR_ARGS(ChasteCuboid<1>, p_cuboid, (lower, upper));

        // 0.3 is the MeinekeDivisionSeparation value 
        c_vector<double, 1> expected_daughter_location_vessel_segment  = old_parent_location_forvesselsegment;
        c_vector<double, 1> expected_parent_location_vessel_segment;
        expected_parent_location_vessel_segment[0] = old_parent_location_forvesselsegment[0] - 0.3*1;  

        // we create the sprouting rule object with a MaxSproutingRateConstant of 1.0
        typedef SproutingRuleWithPdes<1,1> SproutingRuleWithPdes;

        // Create PDE 
        MAKE_PTR_ARGS(VegfEquationPde<1>, p_pde, (cell_population, 1.0, 1.0, 1.0, 0.1, 0.01));

        // for sourceterm = 0.1
        MAKE_PTR_ARGS(VegfBoundaryCondition<1>, p_bc_vessel_segment, (0.1, 0));
        MAKE_PTR_ARGS(MolecularConcentrationsDomainPdeModifier<1>, p_pde_modifier_vessel_segment, (p_pde, p_bc_vessel_segment, false, p_cuboid, 1.0,nullptr,0.0,1.0, 0.1));
        p_pde_modifier_vessel_segment->SetDependentVariableName("vegf_femesh_variable");
        p_pde_modifier_vessel_segment->SetOutputGradient(true);
        p_pde_modifier_vessel_segment->SetupInitialSolutionVector(cell_population);

        MAKE_PTR_ARGS(SproutingRuleWithPdes, p_sprouting_rule_vessel_segment, (0.08, 2.0, p_pde_modifier_vessel_segment, 1.0, 0.3, 0.98, 0.4));

        TS_ASSERT_DELTA(p_sprouting_rule_vessel_segment->GetSproutingProbability(cell_population, p_cell_forvesselsegment), 0.001, 1e-4);

        // case of vessel segment 
        // set the division rule to be the division rule for the cell population
        cell_population.SetCentreBasedDivisionRule(p_sprouting_rule_vessel_segment);
        boost::shared_ptr<AbstractCentreBasedDivisionRule<1,1> > p_sprouting_rule_vessel_segment_test = cell_population.GetCentreBasedDivisionRule();

        // test to check that the division rule returns the correct pair of vectors 
        std::pair<c_vector<double, 1>, c_vector<double, 1> > positions_vessel_segment = p_sprouting_rule_vessel_segment_test->CalculateCellDivisionVector(p_cell_forvesselsegment, cell_population);

        c_vector<double, 1> parent_location_vessel_segment;
        parent_location_vessel_segment = positions_vessel_segment.first;
        TS_ASSERT_DELTA(parent_location_vessel_segment[0], expected_parent_location_vessel_segment[0], 1e-6);

        c_vector<double, 1> daughter_location_vessel_segment;
        daughter_location_vessel_segment = positions_vessel_segment.second;
        TS_ASSERT_DELTA(daughter_location_vessel_segment[0], expected_daughter_location_vessel_segment[0], 1e-6);

        // archiving
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "sprouting_rule.arch";
        {
            // Create PDE and boundary condition objects
            MAKE_PTR_ARGS(VegfEquationPde<1>, p_pde, (cell_population, 1.0, 1.0, 1.0, 0.1, 0.01));
            MAKE_PTR_ARGS(VegfBoundaryCondition<1>, p_bc, (0.1, 8));
            MAKE_PTR_ARGS(MolecularConcentrationsDomainPdeModifier<1>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 1.0,nullptr,0.0,0.1, 0.1));
            SproutingRuleWithPdes* const p_sprouting_rule = new SproutingRuleWithPdes(0.08, 2.0, p_pde_modifier, 1, 0.3, 0.98, 0.4);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_sprouting_rule;

            delete p_sprouting_rule;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            SproutingRuleWithPdes* p_sprouting_rule;
            input_arch >> p_sprouting_rule;

            delete p_sprouting_rule;
        }
    }

    void NoTestSproutingRule2WithPdesIn1D() 
    {
        // creation of the mesh
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0u, false, 2));
        nodes.push_back(new Node<1>(1u, false, 1));
        nodes.push_back(new Node<1>(2u, false, 0));

        NodesOnlyMesh<1> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<1> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        cell_population.Update();

        // Set up cell data on the cell population : initialisation for the division 
        for (typename NodeBasedCellPopulation<1>::Iterator cell_it = cell_population.Begin(); 
                                                cell_it != cell_population.End(); ++cell_it )
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_it);
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            c_vector<double, 1> node_coordinates = cell_population.GetLocationOfCellCentre(p_cell);

            double tortuosity = norm_2(node_coordinates);

            double t = SimulationTime::Instance()->GetTime(); 

            (*cell_it)->GetCellData()->SetItem("node_index", node_index);
            (*cell_it)->GetCellData()->SetItem("daughter_type", 0.0); // 0.0 = no daughter cell, 1.0 = daughter cell
            (*cell_it)->GetCellData()->SetItem("BranchNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingPoint", 0.0); // 0.0 = no, 1.0 = yes : necessary to modify the mutation state after sprouting in this modifier 
            (*cell_it)->GetCellData()->SetItem("DivisionNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingLeader", 0.0); // indice of the branching cell leading the branch of the cell considered 
            (*cell_it)->GetCellData()->SetItem("BranchingSegment", 0.0);
            (*cell_it)->GetCellData()->SetItem("ConsecutiveBranches", 1.0); 
            (*cell_it)->GetCellData()->SetItem("BirthTimeCell", t);
            (*cell_it)->GetCellData()->SetItem("TortuosityParent", tortuosity);
            (*cell_it)->GetCellData()->SetItem("TortuosityDaughter", tortuosity);
            (*cell_it)->GetCellData()->SetItem("LoopNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("Anastomosis", 0);
        }

        // create a vector of the expected coordinates of the parent and daughter cells : in this case, the tip cell divides into a vessel segment
        CellPtr p_cell_forvesseltip = cell_population.GetCellUsingLocationIndex(2);
        c_vector<double,1> old_parent_location_forvesseltip = cell_population.GetLocationOfCellCentre(p_cell_forvesseltip);

        ChastePoint<1> lower(0.0);
        ChastePoint<1> upper(10.0);
        MAKE_PTR_ARGS(ChasteCuboid<1>, p_cuboid, (lower, upper));

        // 0.3 is the MeinekeDivisionSeparation value 
        c_vector<double, 1> expected_daughter_location_tip_cell;
        expected_daughter_location_tip_cell[0] = old_parent_location_forvesseltip[0] - 0.5*0.3*(-1) + 0.5*0.3*1;  
        c_vector<double, 1> expected_parent_location_tip_cell;
        expected_parent_location_tip_cell[0] = old_parent_location_forvesseltip[0] + 0.5*0.3*(-1) + 0.5*0.3*1;  

        // we create the sprouting rule object with a MaxSproutingRateConstant of 1.0
        typedef SproutingRuleWithPdes<1,1> SproutingRuleWithPdes;

        // Create PDE 
        MAKE_PTR_ARGS(VegfEquationPde<1>, p_pde, (cell_population, 1.0, 1.0, 1.0, 0.1, 0.01));

        // for sourceterm = 1.0
        MAKE_PTR_ARGS(VegfBoundaryCondition<1>, p_bc_tip_cell, (1.0, 0));
        MAKE_PTR_ARGS(MolecularConcentrationsDomainPdeModifier<1>, p_pde_modifier_tip_cell, (p_pde, p_bc_tip_cell, false, p_cuboid, 1.0,nullptr,0.0,1.0, 0.1));
        p_pde_modifier_tip_cell->SetDependentVariableName("vegf_femesh_variable");
        p_pde_modifier_tip_cell->SetOutputGradient(true);
        p_pde_modifier_tip_cell->SetupInitialSolutionVector(cell_population);
        MAKE_PTR_ARGS(SproutingRuleWithPdes, p_sprouting_rule_tip_cell, (0.98, 2.0, p_pde_modifier_tip_cell, 1, 0.3, 0.98, 0.4));

        TS_ASSERT_DELTA(p_sprouting_rule_tip_cell->GetSproutingProbability(cell_population, p_cell_forvesseltip), 0.97, 1e-2);

        // case of tip cell (take into account new position of the tip cell) 
        // set the division rule to be the division rule for the cell population
        cell_population.SetCentreBasedDivisionRule(p_sprouting_rule_tip_cell);
        boost::shared_ptr<AbstractCentreBasedDivisionRule<1,1> > p_sprouting_rule_tip_cell_test = cell_population.GetCentreBasedDivisionRule();

        // test to check that the division rule returns the correct pair of vectors 
        std::pair<c_vector<double, 1>, c_vector<double, 1> > positions_tip_cell = p_sprouting_rule_tip_cell_test->CalculateCellDivisionVector(p_cell_forvesseltip, cell_population);

        c_vector<double, 1> parent_location_tip_cell;
        parent_location_tip_cell = positions_tip_cell.first;
        TS_ASSERT_DELTA(norm_2(parent_location_tip_cell), norm_2(expected_parent_location_tip_cell), 1e-0);

        c_vector<double, 1> daughter_location_tip_cell;
        daughter_location_tip_cell = positions_tip_cell.second;
        TS_ASSERT_DELTA(norm_2(daughter_location_tip_cell), norm_2(expected_daughter_location_tip_cell), 1e-0);

        // archiving
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "sprouting_rule.arch";
        {
            // Create PDE and boundary condition objects
            MAKE_PTR_ARGS(VegfEquationPde<1>, p_pde, (cell_population, 1.0, 1.0, 1.0, 0.1, 0.01));
            MAKE_PTR_ARGS(VegfBoundaryCondition<1>, p_bc, (1.0, 8));
            MAKE_PTR_ARGS(MolecularConcentrationsDomainPdeModifier<1>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 1.0,nullptr,0.0,1.0, 0.1));
            SproutingRuleWithPdes* const p_sprouting_rule = new SproutingRuleWithPdes(0.98, 2.0, p_pde_modifier, 1, 0.3, 0.98, 0.4);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_sprouting_rule;

            delete p_sprouting_rule;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            SproutingRuleWithPdes* p_sprouting_rule;
            input_arch >> p_sprouting_rule;

            delete p_sprouting_rule;
        }
    }

    void NoTestSproutingRule1WithPdesIn2D() 
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
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

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

        cell_population.Update();

        // Set up cell data on the cell population : initialisation for the division 
        for (typename NodeBasedCellPopulation<2>::Iterator cell_it = cell_population.Begin(); 
                                                cell_it != cell_population.End(); ++cell_it )
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_it);
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            c_vector<double, 2> node_coordinates = cell_population.GetLocationOfCellCentre(p_cell);

            double tortuosity = norm_2(node_coordinates);

            double t = SimulationTime::Instance()->GetTime(); 

            (*cell_it)->GetCellData()->SetItem("node_index", node_index);
            (*cell_it)->GetCellData()->SetItem("daughter_type", 0.0); // 0.0 = no daughter cell, 1.0 = daughter cell
            (*cell_it)->GetCellData()->SetItem("BranchNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingPoint", 0.0); // 0.0 = no, 1.0 = yes : necessary to modify the mutation state after sprouting in this modifier 
            (*cell_it)->GetCellData()->SetItem("DivisionNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingLeader", 0.0); // indice of the branching cell leading the branch of the cell considered 
            (*cell_it)->GetCellData()->SetItem("BranchingSegment", 0.0);
            (*cell_it)->GetCellData()->SetItem("ConsecutiveBranches", 1.0); 
            (*cell_it)->GetCellData()->SetItem("BirthTimeCell", t);
            (*cell_it)->GetCellData()->SetItem("TortuosityParent", tortuosity);
            (*cell_it)->GetCellData()->SetItem("TortuosityDaughter", tortuosity);
            (*cell_it)->GetCellData()->SetItem("LoopNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("Anastomosis", 0);
        }

        // create a vector of the expected coordinates of the parent and daughter cells : in this case, the tip cell divides into a vessel segment
        CellPtr p_cell_forvesselsegment = cell_population.GetCellUsingLocationIndex(2);
        c_vector<double,2> old_parent_location_forvesselsegment = cell_population.GetLocationOfCellCentre(p_cell_forvesselsegment);

        ChastePoint<2> lower(0.0, 0.0);
        ChastePoint<2> upper(10.0, 10.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // 0.3 is the MeinekeDivisionSeparation value 
        c_vector<double, 2> expected_daughter_location_vessel_segment  = old_parent_location_forvesselsegment;
        c_vector<double, 2> expected_parent_location_vessel_segment;
        expected_parent_location_vessel_segment[0] = old_parent_location_forvesselsegment[0] - 0.3*1;  
        expected_parent_location_vessel_segment[1] = old_parent_location_forvesselsegment[1] - 0.3*0;

        // we create the sprouting rule object with a MaxSproutingRateConstant of 1.0
        typedef SproutingRuleWithPdes<2,2> SproutingRuleWithPdes;

        // Create PDE 
        MAKE_PTR_ARGS(VegfEquationPde<2>, p_pde, (cell_population, 1.0, 1.0, 1.0, 0.1, 0.01));

        // for sourceterm = 0.1
        MAKE_PTR_ARGS(VegfBoundaryCondition<2>, p_bc_vessel_segment, (0.1, 0));
        MAKE_PTR_ARGS(MolecularConcentrationsDomainPdeModifier<2>, p_pde_modifier_vessel_segment, (p_pde, p_bc_vessel_segment, false, p_cuboid, 1.0,nullptr,0.0,1.0, 0.1));
        p_pde_modifier_vessel_segment->SetDependentVariableName("vegf_femesh_variable");
        p_pde_modifier_vessel_segment->SetOutputGradient(true);
        p_pde_modifier_vessel_segment->SetupInitialSolutionVector(cell_population);

        MAKE_PTR_ARGS(SproutingRuleWithPdes, p_sprouting_rule_vessel_segment, (0.08, 2.0, p_pde_modifier_vessel_segment, 1.0, 0.3, 0.98, 0.4));

        TS_ASSERT_DELTA(p_sprouting_rule_vessel_segment->GetSproutingProbability(cell_population, p_cell_forvesselsegment), 0.001, 1e-4);

        // case of vessel segment 
        // set the division rule to be the division rule for the cell population
        cell_population.SetCentreBasedDivisionRule(p_sprouting_rule_vessel_segment);
        boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_sprouting_rule_vessel_segment_test = cell_population.GetCentreBasedDivisionRule();

        // test to check that the division rule returns the correct pair of vectors 
        std::pair<c_vector<double, 2>, c_vector<double, 2> > positions_vessel_segment = p_sprouting_rule_vessel_segment_test->CalculateCellDivisionVector(p_cell_forvesselsegment, cell_population);

        c_vector<double, 2> parent_location_vessel_segment;
        parent_location_vessel_segment = positions_vessel_segment.first;
        TS_ASSERT_DELTA(parent_location_vessel_segment[0], expected_parent_location_vessel_segment[0], 1e-6);
        TS_ASSERT_DELTA(parent_location_vessel_segment[1], expected_parent_location_vessel_segment[1], 1e-6);

        c_vector<double, 2> daughter_location_vessel_segment;
        daughter_location_vessel_segment = positions_vessel_segment.second;
        TS_ASSERT_DELTA(daughter_location_vessel_segment[0], expected_daughter_location_vessel_segment[0], 1e-6);
        TS_ASSERT_DELTA(daughter_location_vessel_segment[1], expected_daughter_location_vessel_segment[1], 1e-6);

        // archiving
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "sprouting_rule.arch";
        {
            // Create PDE and boundary condition objects
            MAKE_PTR_ARGS(VegfEquationPde<2>, p_pde, (cell_population, 1.0, 1.0, 1.0, 0.1, 0.01));
            MAKE_PTR_ARGS(VegfBoundaryCondition<2>, p_bc, (0.1, 8));
            MAKE_PTR_ARGS(MolecularConcentrationsDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 1.0,nullptr,0.0,0.1, 0.1));
            SproutingRuleWithPdes* const p_sprouting_rule = new SproutingRuleWithPdes(0.08, 2.0, p_pde_modifier, 1, 0.3, 0.98, 0.4);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_sprouting_rule;

            delete p_sprouting_rule;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            SproutingRuleWithPdes* p_sprouting_rule;
            input_arch >> p_sprouting_rule;

            delete p_sprouting_rule;
        }
    }

    void NoTestSproutingRule2WithPdesIn2D() 
    {
        // creation of the mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0u, false, 2, 5));
        nodes.push_back(new Node<2>(1u, false, 1, 5));
        nodes.push_back(new Node<2>(2u, false, 0, 5));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

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

        cell_population.Update();

        // Set up cell data on the cell population : initialisation for the division 
        for (typename NodeBasedCellPopulation<2>::Iterator cell_it = cell_population.Begin(); 
                                                cell_it != cell_population.End(); ++cell_it )
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_it);
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            c_vector<double, 2> node_coordinates = cell_population.GetLocationOfCellCentre(p_cell);

            double tortuosity = norm_2(node_coordinates);

            double t = SimulationTime::Instance()->GetTime(); 

            (*cell_it)->GetCellData()->SetItem("node_index", node_index);
            (*cell_it)->GetCellData()->SetItem("daughter_type", 0.0); // 0.0 = no daughter cell, 1.0 = daughter cell
            (*cell_it)->GetCellData()->SetItem("BranchNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingPoint", 0.0); // 0.0 = no, 1.0 = yes : necessary to modify the mutation state after sprouting in this modifier 
            (*cell_it)->GetCellData()->SetItem("DivisionNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingLeader", 0.0); // indice of the branching cell leading the branch of the cell considered 
            (*cell_it)->GetCellData()->SetItem("BranchingSegment", 0.0);
            (*cell_it)->GetCellData()->SetItem("ConsecutiveBranches", 1.0); 
            (*cell_it)->GetCellData()->SetItem("BirthTimeCell", t);
            (*cell_it)->GetCellData()->SetItem("TortuosityParent", tortuosity);
            (*cell_it)->GetCellData()->SetItem("TortuosityDaughter", tortuosity);
            (*cell_it)->GetCellData()->SetItem("LoopNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("Anastomosis", 0);
        }

        // create a vector of the expected coordinates of the parent and daughter cells : in this case, the tip cell divides into a vessel segment
        CellPtr p_cell_forvesseltip = cell_population.GetCellUsingLocationIndex(2);
        c_vector<double,2> old_parent_location_forvesseltip = cell_population.GetLocationOfCellCentre(p_cell_forvesseltip);

        ChastePoint<2> lower(0.0, 0.0);
        ChastePoint<2> upper(10.0, 10.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // 0.3 is the MeinekeDivisionSeparation value 
        c_vector<double, 2> expected_daughter_location_tip_cell;
        expected_daughter_location_tip_cell[0] = old_parent_location_forvesseltip[0] - 0.5*0.3*0 + 0.5*0.3*1;  
        expected_daughter_location_tip_cell[1] = old_parent_location_forvesseltip[1] - 0.5*0.3*1 + 0.5*0.3*0;
        c_vector<double, 2> expected_parent_location_tip_cell;
        expected_parent_location_tip_cell[0] = old_parent_location_forvesseltip[0] + 0.5*0.3*0 + 0.5*0.3*1;  
        expected_parent_location_tip_cell[1] = old_parent_location_forvesseltip[1] + 0.5*0.3*1 + 0.5*0.3*0;

        // we create the sprouting rule object with a MaxSproutingRateConstant of 1.0
        typedef SproutingRuleWithPdes<2,2> SproutingRuleWithPdes;

        // Create PDE 
        MAKE_PTR_ARGS(VegfEquationPde<2>, p_pde, (cell_population, 1.0, 1.0, 1.0, 0.1, 0.01));

        // for sourceterm = 1.0
        MAKE_PTR_ARGS(VegfBoundaryCondition<2>, p_bc_tip_cell, (1.0, 0));
        MAKE_PTR_ARGS(MolecularConcentrationsDomainPdeModifier<2>, p_pde_modifier_tip_cell, (p_pde, p_bc_tip_cell, false, p_cuboid, 1.0,nullptr,0.0,1.0, 0.1));
        p_pde_modifier_tip_cell->SetDependentVariableName("vegf_femesh_variable");
        p_pde_modifier_tip_cell->SetOutputGradient(true);
        p_pde_modifier_tip_cell->SetupInitialSolutionVector(cell_population);
        MAKE_PTR_ARGS(SproutingRuleWithPdes, p_sprouting_rule_tip_cell, (0.98, 2.0, p_pde_modifier_tip_cell, 1, 0.3, 0.98, 0.4));

        TS_ASSERT_DELTA(p_sprouting_rule_tip_cell->GetSproutingProbability(cell_population, p_cell_forvesseltip), 0.97, 1e-2);

        // case of tip cell (take into account new position of the tip cell) 
        // set the division rule to be the division rule for the cell population
        cell_population.SetCentreBasedDivisionRule(p_sprouting_rule_tip_cell);
        boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_sprouting_rule_tip_cell_test = cell_population.GetCentreBasedDivisionRule();

        // test to check that the division rule returns the correct pair of vectors 
        std::pair<c_vector<double, 2>, c_vector<double, 2> > positions_tip_cell = p_sprouting_rule_tip_cell_test->CalculateCellDivisionVector(p_cell_forvesseltip, cell_population);

        c_vector<double, 2> parent_location_tip_cell;
        parent_location_tip_cell = positions_tip_cell.first;
        TS_ASSERT_DELTA(norm_2(parent_location_tip_cell), norm_2(expected_parent_location_tip_cell), 1e-0);

        c_vector<double, 2> daughter_location_tip_cell;
        daughter_location_tip_cell = positions_tip_cell.second;
        TS_ASSERT_DELTA(norm_2(daughter_location_tip_cell), norm_2(expected_daughter_location_tip_cell), 1e-0);

        // archiving
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "sprouting_rule.arch";
        {
            // Create PDE and boundary condition objects
            MAKE_PTR_ARGS(VegfEquationPde<2>, p_pde, (cell_population, 1.0, 1.0, 1.0, 0.1, 0.01));
            MAKE_PTR_ARGS(VegfBoundaryCondition<2>, p_bc, (1.0, 8));
            MAKE_PTR_ARGS(MolecularConcentrationsDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 1.0,nullptr,0.0,1.0, 0.1));
            SproutingRuleWithPdes* const p_sprouting_rule = new SproutingRuleWithPdes(0.98, 2.0, p_pde_modifier, 1, 0.3, 0.98, 0.4);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_sprouting_rule;

            delete p_sprouting_rule;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            SproutingRuleWithPdes* p_sprouting_rule;
            input_arch >> p_sprouting_rule;

            delete p_sprouting_rule;
        }
    }

    void NoTestSproutingRule1WithPdesIn3D() 
    {
        // creation of the mesh
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u, false, 10, 5, 5));
        nodes.push_back(new Node<3>(1u, false, 9, 5, 5));
        nodes.push_back(new Node<3>(2u, false, 8, 5, 5));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        cell_population.Update();

        // Set up cell data on the cell population : initialisation for the division 
        for (typename NodeBasedCellPopulation<3>::Iterator cell_it = cell_population.Begin(); 
                                                cell_it != cell_population.End(); ++cell_it )
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_it);
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            c_vector<double, 3> node_coordinates = cell_population.GetLocationOfCellCentre(p_cell);

            double tortuosity = norm_2(node_coordinates);

            double t = SimulationTime::Instance()->GetTime(); 

            (*cell_it)->GetCellData()->SetItem("node_index", node_index);
            (*cell_it)->GetCellData()->SetItem("daughter_type", 0.0); // 0.0 = no daughter cell, 1.0 = daughter cell
            (*cell_it)->GetCellData()->SetItem("BranchNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingPoint", 0.0); // 0.0 = no, 1.0 = yes : necessary to modify the mutation state after sprouting in this modifier 
            (*cell_it)->GetCellData()->SetItem("DivisionNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingLeader", 0.0); // indice of the branching cell leading the branch of the cell considered 
            (*cell_it)->GetCellData()->SetItem("BranchingSegment", 0.0);
            (*cell_it)->GetCellData()->SetItem("ConsecutiveBranches", 1.0); 
            (*cell_it)->GetCellData()->SetItem("BirthTimeCell", t);
            (*cell_it)->GetCellData()->SetItem("TortuosityParent", tortuosity);
            (*cell_it)->GetCellData()->SetItem("TortuosityDaughter", tortuosity);
            (*cell_it)->GetCellData()->SetItem("LoopNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("Anastomosis", 0);
        }

        // create a vector of the expected coordinates of the parent and daughter cells : in this case, the tip cell divides into a vessel segment
        CellPtr p_cell_forvesselsegment = cell_population.GetCellUsingLocationIndex(2);
        c_vector<double,3> old_parent_location_forvesselsegment = cell_population.GetLocationOfCellCentre(p_cell_forvesselsegment);

        ChastePoint<3> lower(0.0, 0.0, 0.0);
        ChastePoint<3> upper(10.0, 10.0, 10.0);
        MAKE_PTR_ARGS(ChasteCuboid<3>, p_cuboid, (lower, upper));

        // 0.3 is the MeinekeDivisionSeparation value 
        c_vector<double, 3> expected_daughter_location_vessel_segment  = old_parent_location_forvesselsegment;
        c_vector<double, 3> expected_parent_location_vessel_segment;
        expected_parent_location_vessel_segment[0] = old_parent_location_forvesselsegment[0] - 0.3*1;  
        expected_parent_location_vessel_segment[1] = old_parent_location_forvesselsegment[1] - 0.3*0;
        expected_parent_location_vessel_segment[2] = old_parent_location_forvesselsegment[2] - 0.3*0;

        // we create the sprouting rule object with a MaxSproutingRateConstant of 1.0
        typedef SproutingRuleWithPdes<3,3> SproutingRuleWithPdes;

        // Create PDE 
        MAKE_PTR_ARGS(VegfEquationPde<3>, p_pde, (cell_population, 1.0, 1.0, 1.0, 0.1, 0.01));

        // for sourceterm = 0.1
        MAKE_PTR_ARGS(VegfBoundaryCondition<3>, p_bc_vessel_segment, (0.1, 0));
        MAKE_PTR_ARGS(MolecularConcentrationsDomainPdeModifier<3>, p_pde_modifier_vessel_segment, (p_pde, p_bc_vessel_segment, false, p_cuboid, 1.0,nullptr,0.0,1.0, 0.1));
        p_pde_modifier_vessel_segment->SetDependentVariableName("vegf_femesh_variable");
        p_pde_modifier_vessel_segment->SetOutputGradient(true);
        p_pde_modifier_vessel_segment->SetupInitialSolutionVector(cell_population);

        MAKE_PTR_ARGS(SproutingRuleWithPdes, p_sprouting_rule_vessel_segment, (0.08, 2.0, p_pde_modifier_vessel_segment, 1.0, 0.3, 0.98, 0.4));

        TS_ASSERT_DELTA(p_sprouting_rule_vessel_segment->GetSproutingProbability(cell_population, p_cell_forvesselsegment), 0.001, 1e-4);

        // case of vessel segment 
        // set the division rule to be the division rule for the cell population
        cell_population.SetCentreBasedDivisionRule(p_sprouting_rule_vessel_segment);
        boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_sprouting_rule_vessel_segment_test = cell_population.GetCentreBasedDivisionRule();

        // test to check that the division rule returns the correct pair of vectors 
        std::pair<c_vector<double, 3>, c_vector<double, 3> > positions_vessel_segment = p_sprouting_rule_vessel_segment_test->CalculateCellDivisionVector(p_cell_forvesselsegment, cell_population);

        c_vector<double, 3> parent_location_vessel_segment;
        parent_location_vessel_segment = positions_vessel_segment.first;
        TS_ASSERT_DELTA(parent_location_vessel_segment[0], expected_parent_location_vessel_segment[0], 1e-6);
        TS_ASSERT_DELTA(parent_location_vessel_segment[1], expected_parent_location_vessel_segment[1], 1e-6);
        TS_ASSERT_DELTA(parent_location_vessel_segment[2], expected_parent_location_vessel_segment[2], 1e-6);

        c_vector<double, 3> daughter_location_vessel_segment;
        daughter_location_vessel_segment = positions_vessel_segment.second;
        TS_ASSERT_DELTA(daughter_location_vessel_segment[0], expected_daughter_location_vessel_segment[0], 1e-6);
        TS_ASSERT_DELTA(daughter_location_vessel_segment[1], expected_daughter_location_vessel_segment[1], 1e-6);
        TS_ASSERT_DELTA(daughter_location_vessel_segment[2], expected_daughter_location_vessel_segment[2], 1e-6);

        // archiving
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "sprouting_rule.arch";
        {
            // Create PDE and boundary condition objects
            MAKE_PTR_ARGS(VegfEquationPde<3>, p_pde, (cell_population, 1.0, 1.0, 1.0, 0.1, 0.01));
            MAKE_PTR_ARGS(VegfBoundaryCondition<3>, p_bc, (0.1, 8));
            MAKE_PTR_ARGS(MolecularConcentrationsDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 1.0,nullptr,0.0,0.1, 0.1));
            SproutingRuleWithPdes* const p_sprouting_rule = new SproutingRuleWithPdes(0.08, 2.0, p_pde_modifier, 1, 0.3, 0.98, 0.4);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_sprouting_rule;

            delete p_sprouting_rule;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            SproutingRuleWithPdes* p_sprouting_rule;
            input_arch >> p_sprouting_rule;

            delete p_sprouting_rule;
        }
    }

    void NoTestSproutingRule2WithPdesIn3D() 
    {
        // creation of the mesh
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u, false, 2, 5, 5));
        nodes.push_back(new Node<3>(1u, false, 1, 5, 5));
        nodes.push_back(new Node<3>(2u, false, 0, 5, 5));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state); 
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide 

        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        cell_population.Update();

        // Set up cell data on the cell population : initialisation for the division 
        for (typename NodeBasedCellPopulation<3>::Iterator cell_it = cell_population.Begin(); 
                                                cell_it != cell_population.End(); ++cell_it )
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_it);
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(node_index);
            c_vector<double, 3> node_coordinates = cell_population.GetLocationOfCellCentre(p_cell);

            double tortuosity = norm_2(node_coordinates);

            double t = SimulationTime::Instance()->GetTime(); 

            (*cell_it)->GetCellData()->SetItem("node_index", node_index);
            (*cell_it)->GetCellData()->SetItem("daughter_type", 0.0); // 0.0 = no daughter cell, 1.0 = daughter cell
            (*cell_it)->GetCellData()->SetItem("BranchNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingPoint", 0.0); // 0.0 = no, 1.0 = yes : necessary to modify the mutation state after sprouting in this modifier 
            (*cell_it)->GetCellData()->SetItem("DivisionNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingLeader", 0.0); // indice of the branching cell leading the branch of the cell considered 
            (*cell_it)->GetCellData()->SetItem("BranchingSegment", 0.0);
            (*cell_it)->GetCellData()->SetItem("ConsecutiveBranches", 1.0); 
            (*cell_it)->GetCellData()->SetItem("BirthTimeCell", t);
            (*cell_it)->GetCellData()->SetItem("TortuosityParent", tortuosity);
            (*cell_it)->GetCellData()->SetItem("TortuosityDaughter", tortuosity);
            (*cell_it)->GetCellData()->SetItem("LoopNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("Anastomosis", 0);
        }

        // create a vector of the expected coordinates of the parent and daughter cells : in this case, the tip cell divides into a vessel segment
        CellPtr p_cell_forvesseltip = cell_population.GetCellUsingLocationIndex(2);
        c_vector<double,3> old_parent_location_forvesseltip = cell_population.GetLocationOfCellCentre(p_cell_forvesseltip);

        ChastePoint<3> lower(0.0, 0.0, 0.0);
        ChastePoint<3> upper(10.0, 10.0, 10.0);
        MAKE_PTR_ARGS(ChasteCuboid<3>, p_cuboid, (lower, upper));

        // 0.3 is the MeinekeDivisionSeparation value 
        c_vector<double, 3> expected_daughter_location_tip_cell;
        expected_daughter_location_tip_cell[0] = old_parent_location_forvesseltip[0] - 0.5*0.3*0 + 0.5*0.3*1;  
        expected_daughter_location_tip_cell[1] = old_parent_location_forvesseltip[1] - 0.5*0.3*1 + 0.5*0.3*0;
        expected_daughter_location_tip_cell[2] = old_parent_location_forvesseltip[2] - 0.5*0.3*0 + 0.5*0.3*0;
        c_vector<double, 3> expected_parent_location_tip_cell;
        expected_parent_location_tip_cell[0] = old_parent_location_forvesseltip[0] + 0.5*0.3*0 + 0.5*0.3*1;  
        expected_parent_location_tip_cell[1] = old_parent_location_forvesseltip[1] + 0.5*0.3*1 + 0.5*0.3*0;
        expected_parent_location_tip_cell[2] = old_parent_location_forvesseltip[2] + 0.5*0.3*0 + 0.5*0.3*0;

        // we create the sprouting rule object with a MaxSproutingRateConstant of 1.0
        typedef SproutingRuleWithPdes<3,3> SproutingRuleWithPdes;

        // Create PDE 
        MAKE_PTR_ARGS(VegfEquationPde<3>, p_pde, (cell_population, 1.0, 1.0, 1.0, 0.1, 0.01));

        // for sourceterm = 1.0
        MAKE_PTR_ARGS(VegfBoundaryCondition<3>, p_bc_tip_cell, (1.0, 0));
        MAKE_PTR_ARGS(MolecularConcentrationsDomainPdeModifier<3>, p_pde_modifier_tip_cell, (p_pde, p_bc_tip_cell, false, p_cuboid, 1.0,nullptr,0.0,1.0, 0.1));
        p_pde_modifier_tip_cell->SetDependentVariableName("vegf_femesh_variable");
        p_pde_modifier_tip_cell->SetOutputGradient(true);
        p_pde_modifier_tip_cell->SetupInitialSolutionVector(cell_population);
        MAKE_PTR_ARGS(SproutingRuleWithPdes, p_sprouting_rule_tip_cell, (0.98, 2.0, p_pde_modifier_tip_cell, 1, 0.3, 0.98, 0.4));

        TS_ASSERT_DELTA(p_sprouting_rule_tip_cell->GetSproutingProbability(cell_population, p_cell_forvesseltip), 0.97, 1e-2);

        // case of tip cell (take into account new position of the tip cell) 
        // set the division rule to be the division rule for the cell population
        cell_population.SetCentreBasedDivisionRule(p_sprouting_rule_tip_cell);
        boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_sprouting_rule_tip_cell_test = cell_population.GetCentreBasedDivisionRule();

        // test to check that the division rule returns the correct pair of vectors 
        std::pair<c_vector<double, 3>, c_vector<double, 3> > positions_tip_cell = p_sprouting_rule_tip_cell_test->CalculateCellDivisionVector(p_cell_forvesseltip, cell_population);

        c_vector<double, 3> parent_location_tip_cell;
        parent_location_tip_cell = positions_tip_cell.first;
        TS_ASSERT_DELTA(norm_2(parent_location_tip_cell), norm_2(expected_parent_location_tip_cell), 1e-0);

        c_vector<double, 3> daughter_location_tip_cell;
        daughter_location_tip_cell = positions_tip_cell.second;
        TS_ASSERT_DELTA(norm_2(daughter_location_tip_cell), norm_2(expected_daughter_location_tip_cell), 1e-0);

        // archiving
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "sprouting_rule.arch";
        {
            // Create PDE and boundary condition objects
            MAKE_PTR_ARGS(VegfEquationPde<3>, p_pde, (cell_population, 1.0, 1.0, 1.0, 0.1, 0.01));
            MAKE_PTR_ARGS(VegfBoundaryCondition<3>, p_bc, (1.0, 8));
            MAKE_PTR_ARGS(MolecularConcentrationsDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 1.0,nullptr,0.0,1.0, 0.1));
            SproutingRuleWithPdes* const p_sprouting_rule = new SproutingRuleWithPdes(0.98, 2.0, p_pde_modifier, 1, 0.3, 0.98, 0.4);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_sprouting_rule;

            delete p_sprouting_rule;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            SproutingRuleWithPdes* p_sprouting_rule;
            input_arch >> p_sprouting_rule;

            delete p_sprouting_rule;
        }
    }

    void NoTestForcesTimeIn2D()
    {
        ///////////
        // INPUT // 
        ///////////

        // parameters for vegf concentration
        double input_val_vegf_constantbackground = 0.5; //0.1

        // parameters for forces 
        double input_val_S = 1E-5; //1E-5
        double input_val_sigma = 1E-1; //1E-1
        double input_val_chi = 1E-4; //1E-4
        double input_val_omegap = 1E-1; //1E-1
        double input_val_omegaa = 1E-5; //1E-5

        // parameters for Psprout 
        // double input_val_maxsproutingrate = 0.08; 
        // double input_val_cmax = 0.8;
        // double input_val_cmin = 0.3;
        // double input_val_pmax = 0.98;
        // double input_val_pmin = 0.5;

        // parameters for anastomosis 
        // double input_val_anastomosislength = 1.5;
        // double input_val_thresholdlength = 2.0;

        // general parameters (time, random seed, output directory)
        double input_val_time = 1; 
        double input_val_seed = 10; 
        std::string output_directory = "PaperAngiogenesisModel2025/PaperModel2025Analysis2D/CoupledModel2DConstant/Test6/CoupledModel2DConstantSeed10Source5 -output_directory_twolesions PaperAngiogenesisModel2025/PaperModel2025Analysis2D/CoupledModel2DTwoLesions/CoupledModel2DConstantSeed10Sources5OneStep";

        ////////////
        // MESHES // 
        ////////////

        // mesh for pdes 

        // creation of the mesh for ABM
        // read results.viznodes to obtain the coordinates of the centre of the cells 
        std::ifstream infile_nodes;
        std::string line_nodes;
        std::string last_line_nodes;

        infile_nodes.open("/home/chaste/testoutput/PaperAngiogenesisModel2025/PaperModel2025Analysis2D/CoupledModel2DConstant/Test6/CoupledModel2DConstantSeed10Source5/results_from_time_0/results.viznodes");

        if (!infile_nodes.is_open())
        {
            std::cerr << "Error opening file: " << "results.viznodes" << std::endl;
        }

        // Get the last non-empty line
        while (std::getline(infile_nodes, line_nodes))
        {
            if (!line_nodes.empty())
            {
                last_line_nodes = line_nodes;
            }
        }

        if (last_line_nodes.empty())
        {
            EXCEPTION("Input file is empty or improperly formatted.");
        }

        std::stringstream ss_nodes(last_line_nodes);

        double timestep;
        ss_nodes >> timestep;  // Discard the first value (timestep)

        std::vector<Node<2>*> nodes;
        double x, y;
        unsigned index = 0;

        // Parse the remaining coordinates as node positions
        while (index < 10507)
        {
            ss_nodes >> x >> y;

            c_vector<double, 2> location;
            location[0] = x;
            location[1] = y;

            Node<2>* new_node = new Node<2>(index, location, false);
            nodes.push_back(new_node);
            ++index;
        }

        infile_nodes.close();

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        ///////////
        // CELLS // 
        ///////////

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(VesselTipMutationState, p_tip_state);
        MAKE_PTR(VesselSegmentMutationState, p_vessel_state);
        MAKE_PTR(BranchingSegmentMutationState, p_branching_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide, fully constrained 

        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 10507, p_differentiated_type);

        // read results.vizmutationstates to obtain the cell 
        std::ifstream infile_mutations;
        std::string line_mutations;
        std::string last_line_mutations;

        infile_mutations.open("/home/chaste/testoutput/PaperAngiogenesisModel2025/PaperModel2025Analysis2D/CoupledModel2DConstant/Test6/CoupledModel2DConstantSeed10Source5/results_from_time_0/results.vizmutationstates");

        // Get the last non-empty line
        while (std::getline(infile_mutations, line_mutations))
        {
            if (!line_mutations.empty())
            {
                last_line_mutations = line_mutations;
            }
        }

        if (last_line_mutations.empty())
        {
            EXCEPTION("Input file is empty or improperly formatted.");
        }

        std::stringstream ss_mutations(last_line_mutations);

        double timestep_mutations;
        ss_mutations >> timestep_mutations;  // Discard the first value (timestep)

        double mutation;
        unsigned index_mutations = 0;

        // Parse the remaining coordinates as node positions
        while (index_mutations < 10507)
        {
            ss_mutations >> mutation;

            if(mutation == 0){
                cells[index_mutations]->SetMutationState(p_tip_state);
            } else if(mutation == 1) {
                cells[index_mutations]->SetMutationState(p_vessel_state);
            } else {
                cells[index_mutations]->SetMutationState(p_branching_state);
            }

            ++index_mutations;
        }

        infile_mutations.close();
        
        // creation of a population of cells 
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Initialise all node forces to zero 
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        cell_population.Update();

        std::vector<std::string> mDataLabels = std::vector<std::string>{ "fx", "fy", "fz" };
        for ( typename NodeBasedCellPopulation<2>::Iterator cell_it = cell_population.Begin(); 
                                                cell_it != cell_population.End(); ++cell_it )
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_it);

            for (unsigned i = 0; i < 2; i++)
            {
                (*cell_it)->GetCellData()->SetItem(mDataLabels[i], 0.0);
            }

            (*cell_it)->GetCellData()->SetItem("node_index", node_index);
            (*cell_it)->GetCellData()->SetItem("daughter_type", 0.0); // 0.0 = no daughter cell, 1.0 = daughter cell
            (*cell_it)->GetCellData()->SetItem("BranchNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingPoint", 0.0); // 0.0 = no, 1.0 = yes : necessary to modify the mutation state after sprouting in this modifier 
            (*cell_it)->GetCellData()->SetItem("DivisionNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("BranchingLeader", 0.0); // indice of the branching cell leading the branch of the cell considered 
            (*cell_it)->GetCellData()->SetItem("BranchingSegment", 0.0);
            (*cell_it)->GetCellData()->SetItem("ConsecutiveBranches", 1.0); 
            (*cell_it)->GetCellData()->SetItem("BirthTimeCell", 0.0);
            (*cell_it)->GetCellData()->SetItem("TortuosityParent", 0.0);
            (*cell_it)->GetCellData()->SetItem("TortuosityDaughter", 0.0);
            (*cell_it)->GetCellData()->SetItem("LoopNumber", 0.0);
            (*cell_it)->GetCellData()->SetItem("Anastomosis", 0);
        }

        //////////////
        // SETTINGS // 
        //////////////

        // set seed 
        RandomNumberGenerator::Instance()->Reseed(input_val_seed);

        // Set up simulation time for file output
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(input_val_time);

        /////////////////
        // SIMULATION // 
        ////////////////

        // 1) SOLVE PDE 

        // No PDE in the case of constant VEGF in the ECM

        // 2) UPDATING CELL POSITION 

        // Random force (all cells)
        typedef RandomForce<2> RandomForce;
        MAKE_PTR_ARGS(RandomForce, p_random_force, (input_val_sigma));
        simulator.AddForce(p_random_force);

        // Chemotactic force (tip cells only) 
        typedef ChemoForceWithConstantVegf<2> ChemoForceWithConstantVegf;
        MAKE_PTR_ARGS(ChemoForceWithConstantVegf, p_chemo_force, (input_val_chi, input_val_vegf_constantbackground)); 
        simulator.AddForce(p_chemo_force);

        //Persistence force (tip cells only)
        typedef PersistenceForce<2> PersistenceForce;
        MAKE_PTR_ARGS(PersistenceForce, p_persistence_force, (input_val_omegap));
        simulator.AddForce(p_persistence_force);

        // Mechanical force (all cells)
        typedef LinearMechanicalForceModified<2> LinearMechanicalForceModified;
        MAKE_PTR(LinearMechanicalForceModified, p_mechanical_force);
        p_mechanical_force->SetMeinekeSpringStiffness(15.0*input_val_S);
        cell_population.SetMeinekeDivisionSeparation(1.0);
        p_mechanical_force->SetMeinekeDivisionRestingSpringLength(1.0);
        p_mechanical_force->SetMeinekeSpringGrowthDuration(1.0);
        p_mechanical_force->SetCutOffLength(1.5);
        simulator.AddForce(p_mechanical_force);

        // Angular force (vessel segment only)
        typedef AngularForce<2> AngularForce;
        MAKE_PTR_ARGS(AngularForce, p_angular_force, (-input_val_omegaa)); 
        simulator.AddForce(p_angular_force);

        simulator.Solve();

        // Output run time data
        CellBasedEventHandler::Headings();
        CellBasedEventHandler::Report();

        SimulationTime::Destroy();
    }
};

#endif /*TESTFORCESMODEL_HPP_*/