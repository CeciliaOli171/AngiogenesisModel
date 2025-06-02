#ifndef TESTANGIOGENESISMODELWITHVEGFCONCENTRATIONPDE_HPP_
#define TESTANGIOGENESISMODELWITHVEGFCONCENTRATIONPDE_HPP_

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
#include "EllipticBoxDomainPdeModifier.hpp"
#include "ParabolicBoxDomainPdeModifier.hpp"
#include "AveragedSourceEllipticPde.hpp"
#include "AveragedSourceParabolicPde.hpp"
#include "UniformSourceParabolicPde.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "MolecularConcentrationsDomainPdeModifier.hpp"
#include "CellwiseSourceParabolicPde.hpp"
#include "VegfBoundaryCondition.hpp"



class TestAngiogenesisModelWithVegfConcentrationPdes : public AbstractCellBasedTestSuite
{

public:

    void NoTestAngiogenesisModelWithVegfConcentrationPdeIn2D() 
    {
        ///////////
        // INPUT // 
        ///////////

        // Parameters input 
        CommandLineArguments* command_line = CommandLineArguments::Instance();

        // parameters for vegf pde 
        double input_val_vegf_dudtcoeff = command_line->GetDoubleCorrespondingToOption("-vegfdudtcoeff"); //1.0
        double input_val_vegf_diffusioncoeff = command_line->GetDoubleCorrespondingToOption("-vegfdiffusioncoeff"); //0.2
        double input_val_vegf_decaycoeff = command_line->GetDoubleCorrespondingToOption("-vegfdecaycoeff"); //1.0
        double input_val_vegf_creationcoeff = command_line->GetDoubleCorrespondingToOption("-vegfcreationcoeff"); //0.0
        double input_val_vegf_consumptioncoeff = command_line->GetDoubleCorrespondingToOption("-vegfconsumptioncoeff"); //0.0
        double input_val_vegf_initialvalue = command_line->GetDoubleCorrespondingToOption("-vegfinitialvalue"); //0.1
        double input_val_vegf_boundaryvalue = command_line->GetDoubleCorrespondingToOption("-vegfboundaryvalue"); //0.1
        double input_val_vegf_constantbackground = command_line->GetDoubleCorrespondingToOption("-vegfconstantbackground"); //0.1

        // parameters for forces 
        double input_val_sigma = command_line->GetDoubleCorrespondingToOption("-sigma"); //1E-4
        double input_val_chi = command_line->GetDoubleCorrespondingToOption("-chi"); //1E-1
        double input_val_omegap = command_line->GetDoubleCorrespondingToOption("-omegap"); //1E-4
        double input_val_omegaa = command_line->GetDoubleCorrespondingToOption("-omegaa"); //1.0

        // parameters for Psprout 
        double input_val_maxsproutingrate = command_line->GetDoubleCorrespondingToOption("-maxsproutingrate"); 
        int input_psproutfunctiontestnb = command_line->GetIntCorrespondingToOption("-psproutfunctiontestnb");

        // parameters for anastomosis 
        double input_val_anastomosislength = commandline->GetDoubleCorrespondingToOption("-anastomosislength");

        // general parameters (time, random seed, output directory)
        double input_val_time = command_line->GetDoubleCorrespondingToOption("-time"); 
        double input_val_seed = command_line->GetIntCorrespondingToOption("-seed"); 
        std::string output_directory = command_line->GetStringCorrespondingToOption("-output_directory");

        ////////////
        // MESHES // 
        ////////////

        // mesh for pdes 
        double boundary_cuboid_min = 0.0;
        double boundary_cuboid_max = 100.0;

        // creation of the mesh for ABM
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0u, false, boundary_cuboid_max*0.25, boundary_cuboid_max/2));
        nodes.push_back(new Node<2>(1u, false, boundary_cuboid_max*0.25-1, boundary_cuboid_max/2));
        nodes.push_back(new Node<2>(2u, false, boundary_cuboid_max*0.25-2, boundary_cuboid_max/2));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        ///////////
        // CELLS // 
        ///////////

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingCellMutationState, p_branching_state); 
        MAKE_PTR(TipCellMutationState, p_tip_state);
        MAKE_PTR(VesselCellMutationState, p_vessel_state);

        // proliferative states
        MAKE_PTR(StemCellProliferativeType, p_stem_type); // all cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type); // first cell cannot divide, fully constrained 

        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_differentiated_type);

        cells[2]->SetCellProliferativeType(p_stem_type);

        cells[0]->SetMutationState(p_vessel_state);
        cells[1]->SetMutationState(p_vessel_state);
        cells[2]->SetMutationState(p_tip_state);
        
        // creation of a population of cells 
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        cell_population.Update(); // addition of this line compared to the sprouting test

        // Writers 
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<ConsecutiveBranchesWriter>();
        cell_population.AddCellWriter<BranchesNumberWriter>();
        cell_population.AddCellWriter<BirthTimeCellWriter>();
        cell_population.AddCellWriter<TortuosityWriter>();

        // fully constrain the first cell using the boundary condition 
        unsigned node_index_tip_cell = cell_population.GetLocationIndexUsingCell(0);
        std::vector<unsigned> pinned_node_indices;
        pinned_node_indices.push_back(node_index_tip_cell);
        typedef PinnedCellsBoundaryCondition<2,2> PinnedCellsBoundaryCondition;
        MAKE_PTR_ARGS(PinnedCellsBoundaryCondition, p_boundary_condition, (&cell_population, pinned_node_indices));

        //////////////
        // SETTINGS // 
        //////////////

        // set seed 
        RandomNumberGenerator::Instance()->Reseed(input_val_seed);

        // Set up simulation time for file output
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(input_val_time);
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

        /////////////////
        // SIMULATION // 
        ////////////////

        // 1) SOLVE PDE 

        // Create PDE and boundary condition objects
        typedef VegfEquationPde<2> VegfEquationPde; 
        typedef VegfBoundaryCondition<2> VegfBoundaryCondition;
        MAKE_PTR_ARGS(VegfEquationPde, p_vegf_pde, (cell_population, input_val_vegf_dudtcoeff, input_val_vegf_diffusioncoeff, input_val_vegf_decaycoeff, input_val_vegf_creationcoeff, input_val_vegf_consumptioncoeff, input_val_vegf_constantbackground));
        MAKE_PTR_ARGS(VegfBoundaryCondition, p_vegf_bc, (input_val_vegf_boundaryvalue, boundary_cuboid_min));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(boundary_cuboid_min, boundary_cuboid_min);
        ChastePoint<2> upper(boundary_cuboid_max, boundary_cuboid_max);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Initial conditions
        Vec initial_condition = nullptr;

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        typedef MolecularConcentrationsDomainPdeModifier<2> MolecularConcentrationsDomainPdeModifier;
        MAKE_PTR_ARGS(MolecularConcentrationsDomainPdeModifier, p_pde_modifier, (p_vegf_pde, p_vegf_bc, false, p_cuboid, 1.0, initial_condition, boundary_cuboid_min, input_val_vegf_initialvalue, input_val_vegf_constantbackground));
        p_pde_modifier->SetDependentVariableName("vegf_femesh_variable");
        p_pde_modifier->SetOutputGradient(true);

        simulator.AddSimulationModifier(p_pde_modifier);


        // 2) UPDATING CELL POSITION 

        // Random force (all cells)
        typedef RandomForce<2> RandomForce;
        MAKE_PTR_ARGS(RandomForce, p_random_force, (input_val_sigma));
        simulator.AddForce(p_random_force);

        // Chemotactic force (tip cells only) 
        typedef ChemoForceWithPdes<2> ChemoForceWithPdes;
        MAKE_PTR_ARGS(ChemoForceWithPdes, p_chemo_force, (input_val_chi, p_pde_modifier));
        simulator.AddForce(p_chemo_force);

        //Persistence force (tip cells only)
        typedef PersistenceForce<2> PersistenceForce;
        MAKE_PTR_ARGS(PersistenceForce, p_persistence_force, (input_val_omegap));
        simulator.AddForce(p_persistence_force);

        // Mechanical force (all cells)
        typedef LinearMechanicalForceModified<2> LinearMechanicalForceModified;
        MAKE_PTR(LinearMechanicalForceModified, p_mechanical_force);
        p_mechanical_force->SetMeinekeSpringStiffness(15.0);
        cell_population.SetMeinekeDivisionSeparation(1.0);
        p_mechanical_force->SetMeinekeDivisionRestingSpringLength(1.0);
        p_mechanical_force->SetMeinekeSpringGrowthDuration(1.0);
        p_mechanical_force->SetCutOffLength(1.5);
        simulator.AddForce(p_mechanical_force);

        // Angular force (vessel segment only)
        typedef AngularForce<2> AngularForce;
        MAKE_PTR_ARGS(AngularForce, p_angular_force, (-input_val_omegaa)); 
        simulator.AddForce(p_angular_force);

        // 3) DIVISION OF CELLS 

        // Set the division rule for our population to be the random direction division rule
        typedef SproutingRuleWithPdes<2,2> SproutingRuleWithPdes;
        MAKE_PTR_ARGS(SproutingRuleWithPdes, p_division_rule_to_set, (input_val_maxsproutingrate, p_pde_modifier, input_psproutfunctiontestnb));
        
        // Set the division rule for our population to be the new division rule implemented earlier 
        cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);

        // we set for each new daughter cell in the population if it is a tip cell or a vessel segment by using the function DaughterTypeofCell
        typedef DaughterCellModifier<2> DaughterCellModifier;
        MAKE_PTR_ARGS(DaughterCellModifier, p_daughtercell_modifier, (input_val_anastomosislength));
        simulator.AddSimulationModifier(p_daughtercell_modifier);

        typedef DirectionalPersistenceCellModifier<2> DirectionalPersistenceCellModifier;
        MAKE_PTR_ARGS(DirectionalPersistenceCellModifier, p_persistenceforce_modifier, ());
        simulator.AddSimulationModifier(p_persistenceforce_modifier);

        cell_population.Update();

        simulator.Solve();

        // Output run time data
        CellBasedEventHandler::Headings();
        CellBasedEventHandler::Report();

        SimulationTime::Destroy();
    }

    void TestAngiogenesisModelWithVegfConcentrationAnalyticalApproximationOfPdeIn2D() throw(Exception)
    {
        ///////////
        // INPUT // 
        ///////////

        // Parameters input 
        CommandLineArguments* command_line = CommandLineArguments::Instance();

        // parameters for vegf pde 
        //double input_val_vegf_dudtcoeff = command_line->GetDoubleCorrespondingToOption("-vegfdudtcoeff"); //1.0 but not necessary in the case of analytical approximation 
        double input_val_vegf_diffusioncoeff = command_line->GetDoubleCorrespondingToOption("-vegfdiffusioncoeff"); //0.2
        double input_val_vegf_decaycoeff = command_line->GetDoubleCorrespondingToOption("-vegfdecaycoeff"); //1.0
        double input_val_vegf_creationcoeff = command_line->GetDoubleCorrespondingToOption("-vegfcreationcoeff"); //0.0
        double input_val_vegf_consumptioncoeff = command_line->GetDoubleCorrespondingToOption("-vegfconsumptioncoeff"); //0.0
        //double input_val_vegf_initialvalue = command_line->GetDoubleCorrespondingToOption("-vegfinitialvalue"); //0.1 but not necessary in the case of analytical approximation
        double input_val_vegf_boundaryvalue = command_line->GetDoubleCorrespondingToOption("-vegfboundaryvalue"); //0.1
        double input_val_vegf_constantbackground = command_line->GetDoubleCorrespondingToOption("-vegfconstantbackground"); //0.1

        // parameters for forces 
        double input_val_sigma = command_line->GetDoubleCorrespondingToOption("-sigma"); //1E-4
        double input_val_chi = command_line->GetDoubleCorrespondingToOption("-chi"); //1E-1
        double input_val_omegap = command_line->GetDoubleCorrespondingToOption("-omegap"); //1E-4
        double input_val_omegaa = command_line->GetDoubleCorrespondingToOption("-omegaa"); //1.0

        // parameters for Psprout 
        double input_val_maxsproutingrate = command_line->GetDoubleCorrespondingToOption("-maxsproutingrate"); 
        int input_psproutfunctiontestnb = command_line->GetIntCorrespondingToOption("-psproutfunctiontestnb"); // 0 for linear function; 1 for hill function
        double input_val_cmax = commandline->GetDoubleCorrespondingToOption("-cmax");
        double input_val_cmin = commandline->GetDoubleCorrespondingToOption("-cmin");
        double input_val_pmax = commandline->GetDoubleCorrespondingToOption("-pmax");
        double input_val_pmin = commandline->GetDoubleCorrespondingToOption("-pmin");

        // parameters for anastomosis 
        double input_val_anastomosislength = commandline->GetDoubleCorrespondingToOption("-anastomosislength");

        // general parameters (time, random seed, output directory)
        double input_val_time = command_line->GetDoubleCorrespondingToOption("-time"); 
        double input_val_seed = command_line->GetIntCorrespondingToOption("-seed"); 
        std::string output_directory_analyticalapproxvegf = command_line->GetStringCorrespondingToOption("-output_directory_analyticalapproxvegf");

        ////////////
        // MESHES // 
        ////////////

        // mesh for pdes 
        double boundary_cuboid_min = 0.0;
        double boundary_cuboid_max = 220.0;

        // creation of the mesh for ABM
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0u, false, boundary_cuboid_max-20, boundary_cuboid_max/2));
        nodes.push_back(new Node<2>(1u, false, boundary_cuboid_max-21, boundary_cuboid_max/2));
        nodes.push_back(new Node<2>(2u, false, boundary_cuboid_max-22, boundary_cuboid_max/2));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        ///////////
        // CELLS // 
        ///////////

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

        cell_population.Update(); // addition of this line compared to the sprouting test

        // Writers 
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<ConsecutiveBranchesWriter>();
        cell_population.AddCellWriter<BranchesNumberWriter>();
        cell_population.AddCellWriter<BirthTimeCellWriter>();
        cell_population.AddCellWriter<TortuosityWriter>();

        // fully constrain the first cell using the boundary condition 
        unsigned node_index_tip_cell = cell_population.GetLocationIndexUsingCell(0);
        std::vector<unsigned> pinned_node_indices;
        pinned_node_indices.push_back(node_index_tip_cell);
        typedef PinnedCellsBoundaryCondition<2,2> PinnedCellsBoundaryCondition;
        MAKE_PTR_ARGS(PinnedCellsBoundaryCondition, p_boundary_condition, (&cell_population, pinned_node_indices));

        //////////////
        // SETTINGS // 
        //////////////

        // set seed 
        RandomNumberGenerator::Instance()->Reseed(input_val_seed);

        // Set up simulation time for file output
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory_analyticalapproxvegf);
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(input_val_time);
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

        /////////////////
        // SIMULATION // 
        ////////////////

        // 1) SOLVE PDE 

        // We approximate the solution 


        // 2) UPDATING CELL POSITION 

        // Random force (all cells)
        typedef RandomForce<2> RandomForce;
        MAKE_PTR_ARGS(RandomForce, p_random_force, (input_val_sigma));
        simulator.AddForce(p_random_force);

        // Chemotactic force (tip cells only) 
        typedef ChemoForceWithAnalyticalApproximationPde<2> ChemoForceWithAnalyticalApproximationPde;
        MAKE_PTR_ARGS(ChemoForceWithAnalyticalApproximationPde, p_chemo_force, (input_val_chi, input_val_vegf_diffusioncoeff, input_val_vegf_decaycoeff, input_val_vegf_creationcoeff, input_val_vegf_consumptioncoeff, input_val_vegf_boundaryvalue, 
        input_val_vegf_constantbackground, boundary_cuboid_max, input_val_vegf_boundaryvalue));
        simulator.AddForce(p_chemo_force);

        //Persistence force (tip cells only)
        typedef PersistenceForce<2> PersistenceForce;
        MAKE_PTR_ARGS(PersistenceForce, p_persistence_force, (input_val_omegap));
        simulator.AddForce(p_persistence_force);

        // Mechanical force (all cells)
        typedef LinearMechanicalForceModified<2> LinearMechanicalForceModified;
        MAKE_PTR(LinearMechanicalForceModified, p_mechanical_force);
        p_mechanical_force->SetMeinekeSpringStiffness(15.0);
        cell_population.SetMeinekeDivisionSeparation(1.0);
        p_mechanical_force->SetMeinekeDivisionRestingSpringLength(1.0);
        p_mechanical_force->SetMeinekeSpringGrowthDuration(1.0);
        p_mechanical_force->SetCutOffLength(1.5);
        simulator.AddForce(p_mechanical_force);

        // Angular force (vessel segment only)
        typedef AngularForce<2> AngularForce;
        MAKE_PTR_ARGS(AngularForce, p_angular_force, (-input_val_omegaa)); 
        simulator.AddForce(p_angular_force);

        // 3) DIVISION OF CELLS 

        // Set the division rule for our population to be the random direction division rule
        typedef SproutingRuleWithAnalyticalApproximationPde<2,2> SproutingRuleWithAnalyticalApproximationPde;
        MAKE_PTR_ARGS(SproutingRuleWithAnalyticalApproximationPde, p_division_rule_to_set, (input_val_maxsproutingrate, input_val_vegf_diffusioncoeff, input_val_vegf_decaycoeff, input_val_vegf_creationcoeff, input_val_vegf_consumptioncoeff, input_val_vegf_boundaryvalue, input_val_vegf_constantbackground, input_val_cmax, input_val_cmin, input_val_pmax, input_val_pmin, input_psproutfunctiontestnb));

        // Set the division rule for our population to be the new division rule implemented earlier 
        cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);

        // we set for each new daughter cell in the population if it is a tip cell or a vessel segment by using the function DaughterTypeofCell
        typedef DaughterCellModifier<2> DaughterCellModifier;
        MAKE_PTR_ARGS(DaughterCellModifier, p_daughtercell_modifier, (input_val_anastomosislength));
        simulator.AddSimulationModifier(p_daughtercell_modifier);

        typedef DirectionalPersistenceCellModifier<2> DirectionalPersistenceCellModifier;
        MAKE_PTR_ARGS(DirectionalPersistenceCellModifier, p_persistenceforce_modifier, ());
        simulator.AddSimulationModifier(p_persistenceforce_modifier);

        cell_population.Update();

        simulator.Solve();

        // Output run time data
        CellBasedEventHandler::Headings();
        CellBasedEventHandler::Report();

        SimulationTime::Destroy();
    }

    void NoTestAngiogenesisModelWithVegfConcentrationConstantIn2D() 
    {
        ///////////
        // INPUT // 
        ///////////

        // Parameters input 
        CommandLineArguments* command_line = CommandLineArguments::Instance();

        // parameters for vegf concentration
        double input_val_vegf_constantbackground = command_line->GetDoubleCorrespondingToOption("-vegfconstantbackground"); //0.1

        // parameters for forces 
        double input_val_sigma = command_line->GetDoubleCorrespondingToOption("-sigma"); //1E-5
        double input_val_chi = command_line->GetDoubleCorrespondingToOption("-chi"); //1E-2
        double input_val_omegap = command_line->GetDoubleCorrespondingToOption("-omegap"); //1E-3
        double input_val_omegaa = command_line->GetDoubleCorrespondingToOption("-omegaa"); //1.0

        // parameters for Psprout 
        double input_val_maxsproutingrate = command_line->GetDoubleCorrespondingToOption("-maxsproutingrate"); 

        // parameters for anastomosis 
        double input_val_anastomosislength = commandline->GetDoubleCorrespondingToOption("-anastomosislength");

        // general parameters (time, random seed, output directory)
        double input_val_time = command_line->GetDoubleCorrespondingToOption("-time"); 
        double input_val_seed = command_line->GetIntCorrespondingToOption("-seed"); 
        std::string output_directory_vegfconstant = command_line->GetStringCorrespondingToOption("-output_directory_vegfconstant");

        ////////////
        // MESHES // 
        ////////////

        // mesh for pdes 
        //double boundary_cuboid_min = 0.0;
        double boundary_cuboid_max = 300.0;

        // creation of the mesh for ABM
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0u, false, boundary_cuboid_max, boundary_cuboid_max/2));
        nodes.push_back(new Node<2>(1u, false, boundary_cuboid_max-1, boundary_cuboid_max/2));
        nodes.push_back(new Node<2>(2u, false, boundary_cuboid_max-2, boundary_cuboid_max/2));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        ///////////
        // CELLS // 
        ///////////

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingCellMutationState, p_branching_state); 
        MAKE_PTR(TipCellMutationState, p_tip_state);
        MAKE_PTR(VesselCellMutationState, p_vessel_state);

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
        
        // creation of a population of cells 
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        cell_population.Update(); // addition of this line compared to the sprouting test

        // Writers 
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<ConsecutiveBranchesWriter>();
        cell_population.AddCellWriter<BranchesNumberWriter>();
        cell_population.AddCellWriter<BirthTimeCellWriter>();
        cell_population.AddCellWriter<TortuosityWriter>();

        // fully constrain the first cell using the boundary condition 
        unsigned node_index_tip_cell = cell_population.GetLocationIndexUsingCell(0);
        std::vector<unsigned> pinned_node_indices;
        pinned_node_indices.push_back(node_index_tip_cell);
        typedef PinnedCellsBoundaryCondition<2,2> PinnedCellsBoundaryCondition;
        MAKE_PTR_ARGS(PinnedCellsBoundaryCondition, p_boundary_condition, (&cell_population, pinned_node_indices));

        //////////////
        // SETTINGS // 
        //////////////

        // set seed 
        RandomNumberGenerator::Instance()->Reseed(input_val_seed);

        // Set up simulation time for file output
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory_vegfconstant);
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(input_val_time);
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

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
        typedef ChemoForce<2> ChemoForce;
        MAKE_PTR_ARGS(ChemoForce, p_chemo_force, (input_val_chi, input_val_vegf_constantbackground)); 
        simulator.AddForce(p_chemo_force);

        //Persistence force (tip cells only)
        typedef PersistenceForce<2> PersistenceForce;
        MAKE_PTR_ARGS(PersistenceForce, p_persistence_force, (input_val_omegap));
        simulator.AddForce(p_persistence_force);

        // Mechanical force (all cells)
        typedef LinearMechanicalForceModified<2> LinearMechanicalForceModified;
        MAKE_PTR(LinearMechanicalForceModified, p_mechanical_force);
        p_mechanical_force->SetMeinekeSpringStiffness(15.0);
        cell_population.SetMeinekeDivisionSeparation(1.0);
        p_mechanical_force->SetMeinekeDivisionRestingSpringLength(1.0);
        p_mechanical_force->SetMeinekeSpringGrowthDuration(1.0);
        p_mechanical_force->SetCutOffLength(1.5);
        simulator.AddForce(p_mechanical_force);

        // Angular force (vessel segment only)
        typedef AngularForce<2> AngularForce;
        MAKE_PTR_ARGS(AngularForce, p_angular_force, (-input_val_omegaa)); 
        simulator.AddForce(p_angular_force);

        // 3) DIVISION OF CELLS 

        // Set the division rule for our population to be the random direction division rule
        typedef SproutingRule<2,2> SproutingRule;
        MAKE_PTR_ARGS(SproutingRule, p_division_rule_to_set, (input_val_maxsproutingrate*input_val_vegf_constantbackground));
        
        // Set the division rule for our population to be the new division rule implemented earlier 
        cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);

        // we set for each new daughter cell in the population if it is a tip cell or a vessel segment by using the function DaughterTypeofCell
        typedef DaughterCellModifier<2> DaughterCellModifier;
        MAKE_PTR_ARGS(DaughterCellModifier, p_daughtercell_modifier, (input_val_anastomosislength));
        simulator.AddSimulationModifier(p_daughtercell_modifier);

        typedef DirectionalPersistenceCellModifier<2> DirectionalPersistenceCellModifier;
        MAKE_PTR_ARGS(DirectionalPersistenceCellModifier, p_persistenceforce_modifier, ());
        simulator.AddSimulationModifier(p_persistenceforce_modifier);

        cell_population.Update();

        simulator.Solve();

        // Output run time data
        CellBasedEventHandler::Headings();
        CellBasedEventHandler::Report();

        SimulationTime::Destroy();
    }

    void NoTestAngiogenesisModelWithVegfConcentrationPdeIn3D() 
    {
        ///////////
        // INPUT // 
        ///////////

        // Parameters input 
        CommandLineArguments* command_line = CommandLineArguments::Instance();

        // parameters for vegf pde 
        double input_val_vegf_dudtcoeff = command_line->GetDoubleCorrespondingToOption("-vegfdudtcoeff"); //1.0
        double input_val_vegf_diffusioncoeff = command_line->GetDoubleCorrespondingToOption("-vegfdiffusioncoeff"); //0.2
        double input_val_vegf_decaycoeff = command_line->GetDoubleCorrespondingToOption("-vegfdecaycoeff"); //1.0
        double input_val_vegf_creationcoeff = command_line->GetDoubleCorrespondingToOption("-vegfcreationcoeff"); //0.0
        double input_val_vegf_consumptioncoeff = command_line->GetDoubleCorrespondingToOption("-vegfconsumptioncoeff"); //0.0
        double input_val_vegf_initialvalue = command_line->GetDoubleCorrespondingToOption("-vegfinitialvalue"); //0.1
        double input_val_vegf_boundaryvalue = command_line->GetDoubleCorrespondingToOption("-vegfboundaryvalue"); //0.1
        double input_val_vegf_constantbackground = command_line->GetDoubleCorrespondingToOption("-vegfconstantbackground"); //0.1

        // parameters for forces 
        double input_val_sigma = command_line->GetDoubleCorrespondingToOption("-sigma"); //1E-4
        double input_val_chi = command_line->GetDoubleCorrespondingToOption("-chi"); //1E-1
        double input_val_omegap = command_line->GetDoubleCorrespondingToOption("-omegap"); //1E-4
        double input_val_omegaa = command_line->GetDoubleCorrespondingToOption("-omegaa"); //1.0

        // parameters for Psprout 
        double input_val_maxsproutingrate = command_line->GetDoubleCorrespondingToOption("-maxsproutingrate"); 
        int input_psproutfunctiontestnb = command_line->GetIntCorrespondingToOption("-psproutfunctiontestnb");

        // parameters for anastomosis 
        double input_val_anastomosislength = commandline->GetDoubleCorrespondingToOption("-anastomosislength");

        // general parameters (time, random seed, output directory)
        double input_val_time = command_line->GetDoubleCorrespondingToOption("-time"); 
        double input_val_seed = command_line->GetIntCorrespondingToOption("-seed"); 
        std::string output_directory = command_line->GetStringCorrespondingToOption("-output_directory");

        ////////////
        // MESHES // 
        ////////////

        // mesh for pdes 
        double boundary_cuboid_min = 0.0;
        double boundary_cuboid_max = 50.0;

        // creation of the mesh
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u, false, boundary_cuboid_max/2, boundary_cuboid_max/2, boundary_cuboid_max/2));
        nodes.push_back(new Node<3>(1u, false, boundary_cuboid_max/2-1, boundary_cuboid_max/2, boundary_cuboid_max/2));
        nodes.push_back(new Node<3>(2u, false, boundary_cuboid_max/2-2, boundary_cuboid_max/2, boundary_cuboid_max/2));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        ///////////
        // CELLS // 
        ///////////

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingCellMutationState, p_branching_state); 
        MAKE_PTR(TipCellMutationState, p_tip_state);
        MAKE_PTR(VesselCellMutationState, p_vessel_state);

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

        cell_population.Update(); // addition of this line compared to the sprouting test

        // Writers 
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<ConsecutiveBranchesWriter>();
        cell_population.AddCellWriter<BranchesNumberWriter>();
        cell_population.AddCellWriter<BirthTimeCellWriter>();
        cell_population.AddCellWriter<TortuosityWriter>();

        // fully constrain the first cell using the boundary condition 
        unsigned node_index_tip_cell = cell_population.GetLocationIndexUsingCell(0);
        std::vector<unsigned> pinned_node_indices;
        pinned_node_indices.push_back(node_index_tip_cell);
        typedef PinnedCellsBoundaryCondition<3,3> PinnedCellsBoundaryCondition;
        MAKE_PTR_ARGS(PinnedCellsBoundaryCondition, p_boundary_condition, (&cell_population, pinned_node_indices));

        //////////////
        // SETTINGS // 
        //////////////

        // set seed 
        RandomNumberGenerator::Instance()->Reseed(input_val_seed);

        // Set up simulation time for file output
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(input_val_time);
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);


        /////////////////
        // SIMULATION // 
        ////////////////

        // 1) SOLVE PDE 

        // Create PDE and boundary condition objects
        typedef VegfEquationPde<3> VegfEquationPde; 
        typedef VegfBoundaryCondition<3> VegfBoundaryCondition;
        MAKE_PTR_ARGS(VegfEquationPde, p_pde, (cell_population, input_val_vegf_dudtcoeff, input_val_vegf_diffusioncoeff, input_val_vegf_decaycoeff, input_val_vegf_creationcoeff, input_val_vegf_consumptioncoeff, input_val_vegf_constantbackground));
        MAKE_PTR_ARGS(VegfBoundaryCondition, p_bc, (input_val_vegf_boundaryvalue, boundary_cuboid_min));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<3> lower(boundary_cuboid_min, boundary_cuboid_min, boundary_cuboid_min);
        ChastePoint<3> upper(boundary_cuboid_max, boundary_cuboid_max, boundary_cuboid_max);
        MAKE_PTR_ARGS(ChasteCuboid<3>, p_cuboid, (lower, upper));

        // Initial conditions
        Vec initial_condition = nullptr;

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        typedef MolecularConcentrationsDomainPdeModifier<3> MolecularConcentrationsDomainPdeModifier;
        MAKE_PTR_ARGS(MolecularConcentrationsDomainPdeModifier, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 1.0, initial_condition, boundary_cuboid_min, input_val_vegf_initialvalue, input_val_vegf_constantbackground));
        p_pde_modifier->SetDependentVariableName("vegf_femesh_variable");
        p_pde_modifier->SetOutputGradient(true);

        simulator.AddSimulationModifier(p_pde_modifier);


        // 2) UPDATING CELL POSITION 

        // Random force (all cells)
        typedef RandomForce<3> RandomForce;
        MAKE_PTR_ARGS(RandomForce, p_random_force, (input_val_sigma));
        simulator.AddForce(p_random_force);

        // Chemotactic force (tip cells only) 
        typedef ChemoForceWithPdes<3> ChemoForceWithPdes;
        MAKE_PTR_ARGS(ChemoForceWithPdes, p_chemo_force, (input_val_chi, p_pde_modifier));
        simulator.AddForce(p_chemo_force);

        //Persistence force (tip cells only)
        typedef PersistenceForce<3> PersistenceForce;
        MAKE_PTR_ARGS(PersistenceForce, p_persistence_force, (input_val_omegap));
        simulator.AddForce(p_persistence_force);

        // Mechanical force (all cells)
        typedef LinearMechanicalForceModified<3> LinearMechanicalForceModified;
        MAKE_PTR(LinearMechanicalForceModified, p_mechanical_force);
        p_mechanical_force->SetMeinekeSpringStiffness(15.0);
        cell_population.SetMeinekeDivisionSeparation(1.0);
        p_mechanical_force->SetMeinekeDivisionRestingSpringLength(1.0);
        p_mechanical_force->SetMeinekeSpringGrowthDuration(1.0);
        p_mechanical_force->SetCutOffLength(1.5);
        simulator.AddForce(p_mechanical_force);

        // Angular force (vessel segment only)
        typedef AngularForce<3> AngularForce;
        MAKE_PTR_ARGS(AngularForce, p_angular_force, (-input_val_omegaa)); 
        simulator.AddForce(p_angular_force);

        // 3) DIVISION OF CELLS 

        // Set the division rule for our population to be the random direction division rule
        typedef SproutingRuleWithPdes<3,3> SproutingRuleWithPdes;
        MAKE_PTR_ARGS(SproutingRuleWithPdes, p_division_rule_to_set, (input_val_maxsproutingrate, p_pde_modifier, input_psproutfunctiontestnb));

        // Set the division rule for our population to be the new division rule implemented earlier 
        cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);

        // we set for each new daughter cell in the population if it is a tip cell or a vessel segment by using the function DaughterTypeofCell
        typedef DaughterCellModifier<3> DaughterCellModifier;
        MAKE_PTR_ARGS(DaughterCellModifier, p_daughtercell_modifier, (input_val_anastomosislength));
        simulator.AddSimulationModifier(p_daughtercell_modifier);

        typedef DirectionalPersistenceCellModifier<3> DirectionalPersistenceCellModifier;
        MAKE_PTR_ARGS(DirectionalPersistenceCellModifier, p_persistenceforce_modifier, ());
        simulator.AddSimulationModifier(p_persistenceforce_modifier);

        cell_population.Update();

        simulator.Solve();

        // Output run time data
        CellBasedEventHandler::Headings();
        CellBasedEventHandler::Report();

        SimulationTime::Destroy();
    }

    void NoTestAngiogenesisModelWithVegfConcentrationAnalyticalApproximationOfPdeIn3D() 
    {
        ///////////
        // INPUT // 
        ///////////

        // Parameters input 
        CommandLineArguments* command_line = CommandLineArguments::Instance();

        // parameters for vegf pde 
        //double input_val_vegf_dudtcoeff = command_line->GetDoubleCorrespondingToOption("-vegfdudtcoeff"); //1.0 but not necessary for analytical approximation case 
        double input_val_vegf_diffusioncoeff = command_line->GetDoubleCorrespondingToOption("-vegfdiffusioncoeff"); //0.2
        double input_val_vegf_decaycoeff = command_line->GetDoubleCorrespondingToOption("-vegfdecaycoeff"); //1.0
        double input_val_vegf_creationcoeff = command_line->GetDoubleCorrespondingToOption("-vegfcreationcoeff"); //0.0
        double input_val_vegf_consumptioncoeff = command_line->GetDoubleCorrespondingToOption("-vegfconsumptioncoeff"); //0.0
        //double input_val_vegf_initialvalue = command_line->GetDoubleCorrespondingToOption("-vegfinitialvalue"); //0.1 but not necessary for anayltical approximation case 
        double input_val_vegf_boundaryvalue = command_line->GetDoubleCorrespondingToOption("-vegfboundaryvalue"); //0.1
        double input_val_vegf_constantbackground = command_line->GetDoubleCorrespondingToOption("-vegfconstantbackground"); //0.1

        // parameters for forces 
        double input_val_sigma = command_line->GetDoubleCorrespondingToOption("-sigma"); //1E-4
        double input_val_chi = command_line->GetDoubleCorrespondingToOption("-chi"); //1E-1
        double input_val_omegap = command_line->GetDoubleCorrespondingToOption("-omegap"); //1E-4
        double input_val_omegaa = command_line->GetDoubleCorrespondingToOption("-omegaa"); //1.0

        // parameters for Psprout 
        double input_val_maxsproutingrate = command_line->GetDoubleCorrespondingToOption("-maxsproutingrate"); 
        int input_psproutfunctiontestnb = command_line->GetIntCorrespondingToOption("-psproutfunctiontestnb");
        double input_val_cmax = commandline->GetDoubleCorrespondingToOption("-cmax");
        double input_val_cmin = commandline->GetDoubleCorrespondingToOption("-cmin");
        double input_val_pmax = commandline->GetDoubleCorrespondingToOption("-pmax");
        double input_val_pmin = commandline->GetDoubleCorrespondingToOption("-pmin");

        // parameters for anastomosis 
        double input_val_anastomosislength = commandline->GetDoubleCorrespondingToOption("-anastomosislength");

        // general parameters (time, random seed, output directory)
        double input_val_time = command_line->GetDoubleCorrespondingToOption("-time"); 
        double input_val_seed = command_line->GetIntCorrespondingToOption("-seed"); 
        std::string output_directory = command_line->GetStringCorrespondingToOption("-output_directory_analyticalvegf");

        ////////////
        // MESHES // 
        ////////////

        // mesh for pdes 
        //double boundary_cuboid_min = 0.0;
        double boundary_cuboid_max = 100.0;

        // creation of the mesh
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u, false, boundary_cuboid_max*0.25, boundary_cuboid_max/2, boundary_cuboid_max/2));
        nodes.push_back(new Node<3>(1u, false, boundary_cuboid_max*0.25-1, boundary_cuboid_max/2, boundary_cuboid_max/2));
        nodes.push_back(new Node<3>(2u, false, boundary_cuboid_max*0.25-2, boundary_cuboid_max/2, boundary_cuboid_max/2));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        ///////////
        // CELLS // 
        ///////////

        // creation of the cells
        std::vector<CellPtr> cells;

        // mutation states
        MAKE_PTR(BranchingCellMutationState, p_branching_state); 
        MAKE_PTR(TipCellMutationState, p_tip_state);
        MAKE_PTR(VesselCellMutationState, p_vessel_state);

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

        cell_population.Update(); // addition of this line compared to the sprouting test

        // Writers 
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<ConsecutiveBranchesWriter>();
        cell_population.AddCellWriter<BranchesNumberWriter>();
        cell_population.AddCellWriter<BirthTimeCellWriter>();
        cell_population.AddCellWriter<TortuosityWriter>();

        // fully constrain the first cell using the boundary condition 
        unsigned node_index_tip_cell = cell_population.GetLocationIndexUsingCell(0);
        std::vector<unsigned> pinned_node_indices;
        pinned_node_indices.push_back(node_index_tip_cell);
        typedef PinnedCellsBoundaryCondition<3,3> PinnedCellsBoundaryCondition;
        MAKE_PTR_ARGS(PinnedCellsBoundaryCondition, p_boundary_condition, (&cell_population, pinned_node_indices));

        //////////////
        // SETTINGS // 
        //////////////

        // set seed 
        RandomNumberGenerator::Instance()->Reseed(input_val_seed);

        // Set up simulation time for file output
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(input_val_time);
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

        /////////////////
        // SIMULATION // 
        ////////////////

        // 1) SOLVE PDE 

        // We approximate the solution 


        // 2) UPDATING CELL POSITION 

        // Random force (all cells)
        typedef RandomForce<3> RandomForce;
        MAKE_PTR_ARGS(RandomForce, p_random_force, (input_val_sigma));
        simulator.AddForce(p_random_force);

        // Chemotactic force (tip cells only) 
        typedef ChemoForceWithAnalyticalApproximationPde<3> ChemoForceWithAnalyticalApproximationPde;
        MAKE_PTR_ARGS(ChemoForceWithAnalyticalApproximationPde, p_chemo_force, (input_val_chi, input_val_vegf_diffusioncoeff, input_val_vegf_decaycoeff, input_val_vegf_creationcoeff, input_val_vegf_consumptioncoeff, 
        input_val_vegf_boundaryvalue, 
        input_val_vegf_constantbackground, boundary_cuboid_max, input_val_vegf_boundaryvalue));
        simulator.AddForce(p_chemo_force);

        //Persistence force (tip cells only)
        typedef PersistenceForce<3> PersistenceForce;
        MAKE_PTR_ARGS(PersistenceForce, p_persistence_force, (input_val_omegap));
        simulator.AddForce(p_persistence_force);

        // Mechanical force (all cells)
        typedef LinearMechanicalForceModified<3> LinearMechanicalForceModified;
        MAKE_PTR(LinearMechanicalForceModified, p_mechanical_force);
        p_mechanical_force->SetMeinekeSpringStiffness(15.0);
        cell_population.SetMeinekeDivisionSeparation(1.0);
        p_mechanical_force->SetMeinekeDivisionRestingSpringLength(1.0);
        p_mechanical_force->SetMeinekeSpringGrowthDuration(1.0);
        p_mechanical_force->SetCutOffLength(1.5);
        simulator.AddForce(p_mechanical_force);

        // Angular force (vessel segment only)
        typedef AngularForce<3> AngularForce;
        MAKE_PTR_ARGS(AngularForce, p_angular_force, (-input_val_omegaa)); 
        simulator.AddForce(p_angular_force);

        // 3) DIVISION OF CELLS 

        // Set the division rule for our population to be the random direction division rule
        typedef SproutingRuleWithAnalyticalApproximationPde<3,3> SproutingRuleWithAnalyticalApproximationPde;
        MAKE_PTR_ARGS(SproutingRuleWithAnalyticalApproximationPde, p_division_rule_to_set, (input_val_maxsproutingrate, input_val_vegf_diffusioncoeff, input_val_vegf_decaycoeff, input_val_vegf_creationcoeff, input_val_vegf_consumptioncoeff, input_val_vegf_boundaryvalue, input_val_vegf_constantbackground, input_val_cmax, input_val_cmin, input_val_pmax, input_val_pmin, input_psproutfunctiontestnb));

        // Set the division rule for our population to be the new division rule implemented earlier 
        cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);

        // we set for each new daughter cell in the population if it is a tip cell or a vessel segment by using the function DaughterTypeofCell
        typedef DaughterCellModifier<3> DaughterCellModifier;
        MAKE_PTR_ARGS(DaughterCellModifier, p_daughtercell_modifier, (input_val_anastomosislength));
        simulator.AddSimulationModifier(p_daughtercell_modifier);

        typedef DirectionalPersistenceCellModifier<3> DirectionalPersistenceCellModifier;
        MAKE_PTR_ARGS(DirectionalPersistenceCellModifier, p_persistenceforce_modifier, ());
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
        double input_val_sigma = command_line->GetDoubleCorrespondingToOption("-sigma"); //1E-4
        double input_val_chi = command_line->GetDoubleCorrespondingToOption("-chi"); //1E-1
        double input_val_omegap = command_line->GetDoubleCorrespondingToOption("-omegap"); //1E-4
        double input_val_omegaa = command_line->GetDoubleCorrespondingToOption("-omegaa"); //1.0

        // parameters for Psprout 
        double input_val_maxsproutingrate = command_line->GetDoubleCorrespondingToOption("-maxsproutingrate"); 

        // parameters for anastomosis 
        double input_val_anastomosislength = commandline->GetDoubleCorrespondingToOption("-anastomosislength");

        // general parameters (time, random seed, output directory)
        double input_val_time = command_line->GetDoubleCorrespondingToOption("-time"); 
        double input_val_seed = command_line->GetIntCorrespondingToOption("-seed"); 
        std::string output_directory_constantvegf = command_line->GetStringCorrespondingToOption("-output_directory_constantvegf");

        ////////////
        // MESHES // 
        ////////////

        // mesh for pdes 
        //double boundary_cuboid_min = 0.0;
        double boundary_cuboid_max = 100.0;

        // creation of the mesh for ABM
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u, false, boundary_cuboid_max*0.25, boundary_cuboid_max/2, boundary_cuboid_max/2));
        nodes.push_back(new Node<3>(1u, false, boundary_cuboid_max*0.25-1, boundary_cuboid_max/2, boundary_cuboid_max/2));
        nodes.push_back(new Node<3>(2u, false, boundary_cuboid_max*0.25-2, boundary_cuboid_max/2, boundary_cuboid_max/2));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        ///////////
        // CELLS // 
        ///////////

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

        // Writers 
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<ConsecutiveBranchesWriter>();
        cell_population.AddCellWriter<BranchesNumberWriter>();
        cell_population.AddCellWriter<BirthTimeCellWriter>();
        cell_population.AddCellWriter<TortuosityWriter>();

        // fully constrain the first cell using the boundary condition 
        unsigned node_index_tip_cell = cell_population.GetLocationIndexUsingCell(0);
        std::vector<unsigned> pinned_node_indices;
        pinned_node_indices.push_back(node_index_tip_cell);
        typedef PinnedCellsBoundaryCondition<3,3> PinnedCellsBoundaryCondition;
        MAKE_PTR_ARGS(PinnedCellsBoundaryCondition, p_boundary_condition, (&cell_population, pinned_node_indices));

        //////////////
        // SETTINGS // 
        //////////////

        // set seed 
        RandomNumberGenerator::Instance()->Reseed(input_val_seed);

        // Set up simulation time for file output
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory_constantvegf);
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(input_val_time);
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

        /////////////////
        // SIMULATION // 
        ////////////////

        // 1) SOLVE PDE 

        // No PDE in the case of constant VEGF in the ECM 

        // 2) UPDATING CELL POSITION 

        // Random force (all cells)
        typedef RandomForce<3> RandomForce;
        MAKE_PTR_ARGS(RandomForce, p_random_force, (input_val_sigma));
        simulator.AddForce(p_random_force);

        // Chemotactic force (tip cells only) 
        typedef ChemoForce<3> ChemoForce;
        MAKE_PTR_ARGS(ChemoForce, p_chemo_force, (input_val_chi, input_val_vegf_constantbackground));
        simulator.AddForce(p_chemo_force);

        //Persistence force (tip cells only)
        typedef PersistenceForce<3> PersistenceForce;
        MAKE_PTR_ARGS(PersistenceForce, p_persistence_force, (input_val_omegap));
        simulator.AddForce(p_persistence_force);

        // Mechanical force (all cells)
        typedef LinearMechanicalForceModified<3> LinearMechanicalForceModified;
        MAKE_PTR(LinearMechanicalForceModified, p_mechanical_force);
        p_mechanical_force->SetMeinekeSpringStiffness(15.0);
        cell_population.SetMeinekeDivisionSeparation(1.0);
        p_mechanical_force->SetMeinekeDivisionRestingSpringLength(1.0);
        p_mechanical_force->SetMeinekeSpringGrowthDuration(1.0);
        p_mechanical_force->SetCutOffLength(1.5);
        simulator.AddForce(p_mechanical_force);

        // Angular force (vessel segment only)
        typedef AngularForce<3> AngularForce;
        MAKE_PTR_ARGS(AngularForce, p_angular_force, (-input_val_omegaa)); 
        simulator.AddForce(p_angular_force);

        // 3) DIVISION OF CELLS 

        // Set the division rule for our population to be the random direction division rule
        typedef SproutingRule<3,3> SproutingRule;
        MAKE_PTR_ARGS(SproutingRule, p_division_rule_to_set, (input_val_maxsproutingrate*input_val_vegf_constantbackground));

        // Set the division rule for our population to be the new division rule implemented earlier 
        cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);

        // we set for each new daughter cell in the population if it is a tip cell or a vessel segment by using the function DaughterTypeofCell
        typedef DaughterCellModifier<3> DaughterCellModifier;
        MAKE_PTR_ARGS(DaughterCellModifier, p_daughtercell_modifier, (input_val_anastomosislength));
        simulator.AddSimulationModifier(p_daughtercell_modifier);

        typedef DirectionalPersistenceCellModifier<3> DirectionalPersistenceCellModifier;
        MAKE_PTR_ARGS(DirectionalPersistenceCellModifier, p_persistenceforce_modifier, ());
        simulator.AddSimulationModifier(p_persistenceforce_modifier);

        cell_population.Update();

        simulator.Solve();

        // Output run time data
        CellBasedEventHandler::Headings();
        CellBasedEventHandler::Report();

        SimulationTime::Destroy();
    }

};

#endif /*TESTANGIOGENESISMODELWITHVEGFCONCENTRATIONPDE_HPP_*/