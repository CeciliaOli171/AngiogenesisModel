#include "MolecularConcentrationsBoxDomainPdeModifier.hpp"
#include "ParabolicBoxDomainPdeModifier.hpp"
#include "SimpleLinearParabolicSolver.hpp"
#include "VegfEquationPde.hpp"
#include "VegfBoundaryCondition.hpp"
#include "VesselTipMutationState.hpp"
#include "ReplicatableVector.hpp"
#include "AbstractCellPopulation.hpp"
#include "CellLabel.hpp"
#include "CellData.hpp"

#include "Debug.hpp"

template<unsigned DIM>
MolecularConcentrationsBoxDomainPdeModifier<DIM>::MolecularConcentrationsBoxDomainPdeModifier(boost::shared_ptr<AbstractLinearPde<DIM,DIM> > pPde, boost::shared_ptr<AbstractBoundaryCondition<DIM> > pBoundaryCondition, bool isNeumannBoundaryCondition, boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid, double stepSize, Vec solution, double boundaryCuboidMax, double initialValue, double constantBackground)
    : ParabolicBoxDomainPdeModifier<DIM>(pPde, pBoundaryCondition, isNeumannBoundaryCondition, pMeshCuboid, stepSize, solution), mSolutionVegf(solution), mBoundaryCuboidMax(boundaryCuboidMax), mInitialValue(initialValue), mConstantBackground(constantBackground)
{
}

template<unsigned DIM>
MolecularConcentrationsBoxDomainPdeModifier<DIM>::~MolecularConcentrationsBoxDomainPdeModifier()
{
}

template<unsigned DIM>
void MolecularConcentrationsBoxDomainPdeModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Set up boundary conditions
    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc = ConstructBoundaryConditionsContainer(rCellPopulation);

    this->UpdateCellPdeElementMap(rCellPopulation);

    // When using a PDE mesh which doesn't coincide with the cells, we must set up the source terms before solving the PDE.
    // Pass in already updated CellPdeElementMap to speed up finding cells.
    this->SetUpSourceTermsForAveragedSourcePde(this->mpFeMesh, &this->mCellPdeElementMap);
    SimpleLinearParabolicSolver<DIM,DIM> solver(this->mpFeMesh, boost::static_pointer_cast<AbstractLinearParabolicPde<DIM,DIM> >(this->GetPde()).get(), p_bcc.get());

    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double current_time = p_simulation_time->GetTime();
    double dt = p_simulation_time->GetTimeStep();
    solver.SetTimes(current_time,current_time + dt);
    solver.SetTimeStep(dt);

    // Use previous solution as the initial condition
    Vec previous_solution = this->mSolution;
    solver.SetInitialCondition(previous_solution);

    // Note that the linear solver creates a vector, so we have to keep a handle on the old one
    // in order to destroy it
    this->mSolution = solver.Solve();
    PetscTools::Destroy(previous_solution);
    this->UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void MolecularConcentrationsBoxDomainPdeModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    AbstractBoxDomainPdeModifier<DIM>::SetupSolve(rCellPopulation,outputDirectory);
    AbstractBoxDomainPdeModifier<DIM>::UpdateCellData(rCellPopulation);

    // Copy the cell data to mSolution (this is the initial condition)
    ParabolicBoxDomainPdeModifier<DIM>::SetupInitialSolutionVector(rCellPopulation);

    // Output the initial conditions on FeMesh
    this->UpdateAtEndOfOutputTimeStep(rCellPopulation);
}

template<unsigned DIM>
std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > MolecularConcentrationsBoxDomainPdeModifier<DIM>::ConstructBoundaryConditionsContainer(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc(new BoundaryConditionsContainer<DIM,DIM,1>(false));

    if (!this->mSetBcsOnBoxBoundary)
    {
        EXCEPTION("Boundary conditions cannot yet be set on the cell population boundary for a MolecularConcentrationsBoxDomainPdeModifier");
    }
    else // Impose a no-flux Neumann boundary conditions
    {
        for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator elem_iter = this->mpFeMesh->GetBoundaryElementIteratorBegin(); elem_iter != this->mpFeMesh->GetBoundaryElementIteratorEnd(); ++elem_iter)
        {
            p_bcc->AddNeumannBoundaryCondition(*elem_iter, this->mpBoundaryCondition.get());
        }
    }

    return p_bcc;
}


// Explicit instantiation
template class MolecularConcentrationsBoxDomainPdeModifier<1>;
template class MolecularConcentrationsBoxDomainPdeModifier<2>;
template class MolecularConcentrationsBoxDomainPdeModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MolecularConcentrationsBoxDomainPdeModifier)