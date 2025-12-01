#include "MolecularConcentrationsGrowingDomainPdeModifier.hpp"
#include "ParabolicGrowingDomainPdeModifier.hpp"
#include "AveragedSourceParabolicPde.hpp"
#include "CellBasedParabolicPdeSolver.hpp"
#include "VegfEquationCellPde.hpp"
#include "VegfBoundaryCondition.hpp"
#include "VesselTipMutationState.hpp"

#include "Debug.hpp"

template<unsigned DIM>
MolecularConcentrationsGrowingDomainPdeModifier<DIM>::MolecularConcentrationsGrowingDomainPdeModifier(boost::shared_ptr<AbstractLinearPde<DIM,DIM> > pPde, boost::shared_ptr<AbstractBoundaryCondition<DIM> > pBoundaryCondition, bool isNeumannBoundaryCondition, Vec solution, double constantBackground): ParabolicGrowingDomainPdeModifier<DIM>(pPde, pBoundaryCondition, isNeumannBoundaryCondition, solution), mConstantBackground(constantBackground)
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
    UpdateSolutionVectorVEGF(rCellPopulation);

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
    UpdateSolutionVectorVEGF(rCellPopulation);

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

template<unsigned DIM>
void MolecularConcentrationsGrowingDomainPdeModifier<DIM>::UpdateSolutionVectorVEGF(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Clear (if it's not the first time) and resize the solution vector
    if (this->mSolution)
    {
        PetscTools::Destroy(this->mSolution);
    }
    this->mSolution = PetscTools::CreateAndSetVec(this->mpFeMesh->GetNumNodes(), mConstantBackground);

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
void MolecularConcentrationsGrowingDomainPdeModifier<DIM>::SetInitialCellDataVEGF(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin(); cell_iter != rCellPopulation.End(); ++cell_iter)
    {
        cell_iter->GetCellData()->SetItem(this->mDependentVariableName, mConstantBackground);
    }
}

// template<unsigned DIM>
// void MolecularConcentrationsGrowingDomainPdeModifier<DIM>::GenerateFeMeshAngiogenesis(AbstractCellPopulation<DIM,DIM>& rCellPopulation){
//     // initialisation of new cell population 
//     AbstractCellPopulation<DIM,DIM>& rCellPopulationAngiogenesis;

//     // finding the nodes with min and max coordinates : only considerate the new cells and store the old mesh


//     this->GenerateFeMesh(rCellPopulationAngiogenesis);
// }

// Explicit instantiation
template class MolecularConcentrationsGrowingDomainPdeModifier<1>;
template class MolecularConcentrationsGrowingDomainPdeModifier<2>;
template class MolecularConcentrationsGrowingDomainPdeModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MolecularConcentrationsGrowingDomainPdeModifier)