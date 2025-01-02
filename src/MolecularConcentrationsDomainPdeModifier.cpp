#include "MolecularConcentrationsDomainPdeModifier.hpp"
#include "CellBasedParabolicPdeSolver.hpp"
#include "AveragedSourceParabolicPde.hpp"
#include "VegfEquationPde.hpp"

template<unsigned DIM>
MolecularConcentrationsDomainPdeModifier<DIM>::MolecularConcentrationsDomainPdeModifier(boost::shared_ptr<AbstractLinearPde<DIM,DIM> > pPde,
                                                                          boost::shared_ptr<AbstractBoundaryCondition<DIM> > pBoundaryCondition,
                                                                          bool isNeumannBoundaryCondition,
                                                                          Vec solution)
    : AbstractGrowingDomainPdeModifier<DIM>(pPde,
                                            pBoundaryCondition,
                                            isNeumannBoundaryCondition,
                                            solution)
{
}

template<unsigned DIM>
MolecularConcentrationsDomainPdeModifier<DIM>::~MolecularConcentrationsDomainPdeModifier()
{
}

template<unsigned DIM>
void MolecularConcentrationsDomainPdeModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    this->GenerateFeMesh(rCellPopulation);

    // Set up boundary conditions
    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc = ConstructBoundaryConditionsContainer();

    // Construct the solution vector from cell data (takes care of cells dividing);
    UpdateSolutionVector(rCellPopulation);

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
void MolecularConcentrationsDomainPdeModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    AbstractGrowingDomainPdeModifier<DIM>::SetupSolve(rCellPopulation, outputDirectory);

    if (boost::dynamic_pointer_cast<VegfEquationPde<DIM> >(this->mpPde))
    {
        EXCEPTION("MolecularConcentrationsDomainPdeModifier cannot be used with an AveragedSourceParabolicPde. Use a ParabolicBoxDomainPdeModifier instead.");
    }

    // Setup a finite element mesh on which to save the initial condition
    this->GenerateFeMesh(rCellPopulation);

    // Copy the cell data to mSolution (this is the initial condition)
    UpdateSolutionVector(rCellPopulation);

    // Output the initial conditions on FeMesh
    this->UpdateAtEndOfOutputTimeStep(rCellPopulation);
}

template<unsigned DIM>
std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > MolecularConcentrationsDomainPdeModifier<DIM>::ConstructBoundaryConditionsContainer()
{
    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc(new BoundaryConditionsContainer<DIM,DIM,1>(false));

    if (this->IsNeumannBoundaryCondition())
    {
        // Impose any Neumann boundary conditions
        for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator elem_iter = this->mpFeMesh->GetBoundaryElementIteratorBegin();
             elem_iter != this->mpFeMesh->GetBoundaryElementIteratorEnd();
             ++elem_iter)
        {
            p_bcc->AddNeumannBoundaryCondition(*elem_iter, this->mpBoundaryCondition.get());
        }
    }
    else
    {
        // Impose any Dirichlet boundary conditions
        for (typename TetrahedralMesh<DIM,DIM>::BoundaryNodeIterator node_iter = this->mpFeMesh->GetBoundaryNodeIteratorBegin();
             node_iter != this->mpFeMesh->GetBoundaryNodeIteratorEnd();
             ++node_iter)
        {
            p_bcc->AddDirichletBoundaryCondition(*node_iter, this->mpBoundaryCondition.get());
        }
    }

    return p_bcc;
}

template<unsigned DIM>
void MolecularConcentrationsDomainPdeModifier<DIM>::UpdateSolutionVector(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Clear (if it's not the first time) and resize the solution vector
    if (this->mSolution)
    {
        PetscTools::Destroy(this->mSolution);
    }
    this->mSolution = PetscTools::CreateAndSetVec(this->mpFeMesh->GetNumNodes(), 0.0);

    std::string& variable_name = this->mDependentVariableName;

    for (typename TetrahedralMesh<DIM,DIM>::NodeIterator node_iter = this->mpFeMesh->GetNodeIteratorBegin();
         node_iter != this->mpFeMesh->GetNodeIteratorEnd();
         ++node_iter)
    {
        // Loop over nodes of the finite element mesh and get appropriate solution values from CellData
        for (typename TetrahedralMesh<DIM,DIM>::NodeIterator node_iter = this->mpFeMesh->GetNodeIteratorBegin();
             node_iter != this->mpFeMesh->GetNodeIteratorEnd();
             ++node_iter)
        {
            unsigned node_index = node_iter->GetIndex();
            bool dirichlet_bc_applies = (node_iter->IsBoundaryNode()) && (!(this->IsNeumannBoundaryCondition()));
            double boundary_value = this->GetBoundaryCondition()->GetValue(node_iter->rGetLocation());

            double solution_at_node = rCellPopulation.GetCellDataItemAtPdeNode(node_index, variable_name, dirichlet_bc_applies, boundary_value);

            PetscVecTools::SetElement(this->mSolution, node_index, solution_at_node);
        }
    }
}

template<unsigned DIM>
void MolecularConcentrationsDomainPdeModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractGrowingDomainPdeModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class MolecularConcentrationsDomainPdeModifier<1>;
template class MolecularConcentrationsDomainPdeModifier<2>;
template class MolecularConcentrationsDomainPdeModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MolecularConcentrationsDomainPdeModifier)
