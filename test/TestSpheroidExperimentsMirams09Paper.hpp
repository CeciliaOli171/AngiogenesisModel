#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include <iomanip>
#include <boost/foreach.hpp>
#include "OffLatticeSimulation.hpp"
#include "CellBasedEventHandler.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "StochasticOxygenBasedCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "EllipticGrowingDomainPdeModifier.hpp"
#include "ChastePoint.hpp"
#include "SmartPointers.hpp"
#include "ApoptoticCellKiller.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestSpheroidExperimentsMirams98Paper : public AbstractCellBasedTestSuite
{
private:

    void setUp()
    {
        AbstractCellBasedTestSuite::setUp();
        CellBasedEventHandler::Reset();
    }
    void tearDown()
    {
        AbstractCellBasedTestSuite::tearDown();

        CellBasedEventHandler::Headings();
        CellBasedEventHandler::Report();
    }

public:

    void TestMeshBasedSpheroidWithPde() throw(Exception)
    {
        // mesh for the angiogenesis model 
        // std::vector<Node<3>*> nodes;
        // nodes.push_back(new Node<3>(0u, false, 0.0, 0.0, 0.0));
        // nodes.push_back(new Node<3>(1u, false, 1.0, 0.0, 0.0));
        // nodes.push_back(new Node<3>(2u, false, 2.0, 0.0, 0.0));
        // NodesOnlyMesh<3> mesh;
        // mesh.ConstructNodesWithoutMesh(nodes, 1.5); // cut-off length for connectivity of the nodes (=3*Rc=15 for Perfhal's model)

        // original code from Mirams paper 
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true,  1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, true,  1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(3, true,  0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(4, false, 0.5, 0.5, 0.5));
        MutableMesh<3,3> mesh(nodes);

        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        // We use the stem cell G1 duration, so make these 'stem' cells
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<nodes.size(); i++)
        {
            StochasticOxygenBasedCellCycleModel* p_model = new StochasticOxygenBasedCellCycleModel();
            p_model->SetDimension(3);
            p_model->SetStemCellG1Duration(4.0);
            p_model->SetHypoxicConcentration(0.5);
            p_model->SetQuiescentConcentration(0.3);
            p_model->SetCriticalHypoxicDuration(5);
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*
                                (  p_model->GetStemCellG1Duration()
                                 + p_model->GetSG2MDuration() );

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetBirthTime(birth_time);
            p_cell->SetCellProliferativeType(p_stem_type);
            cells.push_back(p_cell);
        }

        MeshBasedCellPopulation<3> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(DBL_MAX);
        cell_population.SetWriteVtkAsPoints(true);

        cell_population.SetDataOnAllCells("oxygen", 1.0);

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetEndTime(100); // hours

        simulator.SetSamplingTimestepMultiple(60);
        simulator.SetOutputDirectory("Plos2013_MeshBasedSpheroidWithPde");

        // Sort out the PDEs
        MAKE_PTR_ARGS(CellwiseSourceEllipticPde<3>, p_pde, (cell_population, -1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("oxygen");

        simulator.AddSimulationModifier(p_pde_modifier);

        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
        p_force->SetMeinekeSpringStiffness(30.0); // default is 15.0;
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(ApoptoticCellKiller<3>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        simulator.Solve();

        CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(&simulator);
    }

    void DNRTestLongerMeshBasedSpheroidWithPde() throw(Exception)
    {
        FileFinder test_data_directory("Plos2013_MeshBasedSpheroidWithPde/archive",
                                       RelativeTo::ChasteTestOutput);

        OutputFileHandler archive_handler("Plos2013_LongerMeshBasedSpheroidWithPde/archive");

        // Following is done in two lines to avoid a bug in Intel compiler v12.0!
        std::vector<FileFinder> temp_files = test_data_directory.FindMatches("*");
        BOOST_FOREACH(FileFinder temp_file, temp_files)
        {
            archive_handler.CopyFileTo(temp_file);
        }

        OffLatticeSimulation<3>* p_simulator
            = CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Load("Plos2013_LongerMeshBasedSpheroidWithPde", 100);

        p_simulator->SetEndTime(150);
        p_simulator->SetOutputDirectory("Plos2013_LongerMeshBasedSpheroidWithPde");

        p_simulator->Solve();

        CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(p_simulator);

    }
};