#include "MolecularConcentrationsGrowingDomainPdeModifier.hpp"
#include "ParabolicGrowingDomainPdeModifier.hpp"
#include "AveragedSourceParabolicPde.hpp"
#include "CellBasedParabolicPdeSolver.hpp"
#include "VegfEquationCellPde.hpp"
#include "VegfBoundaryCondition.hpp"
#include "VesselTipMutationState.hpp"

#include "Debug.hpp"

template<unsigned DIM>
MolecularConcentrationsGrowingDomainPdeModifier<DIM>::MolecularConcentrationsGrowingDomainPdeModifier(boost::shared_ptr<AbstractLinearPde<DIM,DIM> > pPde, boost::shared_ptr<AbstractBoundaryCondition<DIM> > pBoundaryCondition, bool isNeumannBoundaryCondition, Vec solution, double initialValue, double constantBackground): ParabolicGrowingDomainPdeModifier<DIM>(pPde, pBoundaryCondition, isNeumannBoundaryCondition, solution), mInitialValue(initialValue), mConstantBackground(constantBackground)
{
}

template<unsigned DIM>
MolecularConcentrationsGrowingDomainPdeModifier<DIM>::~MolecularConcentrationsGrowingDomainPdeModifier()
{
}

template<unsigned DIM>
void MolecularConcentrationsGrowingDomainPdeModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // generate mesh at each time step to adapt to the cell populaton 
    this->GenerateFeMesh(rCellPopulation);

    // Set up boundary conditions
    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc = ConstructBoundaryConditionsContainer();

    // Construct the solution vector from cell data (takes care of cells dividing);
    ParabolicGrowingDomainPdeModifier<DIM>::UpdateSolutionVector(rCellPopulation);

    // Use CellBasedParabolicPdeSolver as cell wise PDE
    CellBasedParabolicPdeSolver<DIM> solver(this->mpFeMesh, boost::static_pointer_cast<AbstractLinearParabolicPde<DIM,DIM> >(this->mpPde).get(), p_bcc.get());

    ///\todo Investigate more than one PDE time step per spatial step
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
void MolecularConcentrationsGrowingDomainPdeModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    AbstractGrowingDomainPdeModifier<DIM>::SetupSolve(rCellPopulation, outputDirectory);

    if (boost::dynamic_pointer_cast<AveragedSourceParabolicPde<DIM> >(this->mpPde))
    {
        EXCEPTION("MolecularConcentrationsGrowingDomainPdeModifier cannot be used with an AveragedSourceParabolicPde. Use a ParabolicBoxDomainPdeModifier instead.");
    }

    // Setup a finite element mesh on which to save the initial condition
    this->GenerateFeMesh(rCellPopulation);

    // Copy the cell data to mSolution (this is the initial condition)
    ParabolicGrowingDomainPdeModifier<DIM>::UpdateSolutionVector(rCellPopulation);

    // Output the initial conditions on FeMesh
    this->UpdateAtEndOfOutputTimeStep(rCellPopulation);
}

template<unsigned DIM>
std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > MolecularConcentrationsGrowingDomainPdeModifier<DIM>::ConstructBoundaryConditionsContainer()
{
    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc(new BoundaryConditionsContainer<DIM,DIM,1>(false));

    // Impose no-flux Neumann boundary conditions
    for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator elem_iter = this->mpFeMesh->GetBoundaryElementIteratorBegin();elem_iter != this->mpFeMesh->GetBoundaryElementIteratorEnd(); ++elem_iter)
    {
        p_bcc->AddNeumannBoundaryCondition(*elem_iter, this->mpBoundaryCondition.get());
    }

    return p_bcc;
}

// Explicit instantiation
template class MolecularConcentrationsGrowingDomainPdeModifier<1>;
template class MolecularConcentrationsGrowingDomainPdeModifier<2>;
template class MolecularConcentrationsGrowingDomainPdeModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MolecularConcentrationsGrowingDomainPdeModifier)