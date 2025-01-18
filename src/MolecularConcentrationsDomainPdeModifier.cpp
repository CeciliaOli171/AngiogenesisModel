#include "MolecularConcentrationsDomainPdeModifier.hpp"
#include "SimpleLinearParabolicSolver.hpp"
#include "TipCellMutationState.hpp"

#include "Debug.hpp"

template<unsigned DIM>
MolecularConcentrationsDomainPdeModifier<DIM>::MolecularConcentrationsDomainPdeModifier(boost::shared_ptr<AbstractLinearPde<DIM,DIM> > pPde, boost::shared_ptr<AbstractBoundaryCondition<DIM> > pBoundaryCondition, bool isNeumannBoundaryCondition, boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid, double stepSize, Vec solution, double boundaryCuboidMax, double initialValue)
    : AbstractBoxDomainPdeModifier<DIM>(pPde, pBoundaryCondition, isNeumannBoundaryCondition, pMeshCuboid, stepSize, solution), mBoundaryCuboidMax(boundaryCuboidMax), mInitialValue(initialValue)
{
}

template<unsigned DIM>
MolecularConcentrationsDomainPdeModifier<DIM>::~MolecularConcentrationsDomainPdeModifier()
{
}

template<unsigned DIM>
void MolecularConcentrationsDomainPdeModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Set up boundary conditions
    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc = ConstructBoundaryConditionsContainer(rCellPopulation);

    this->UpdateCellPdeElementMap(rCellPopulation);

    // When using a PDE mesh which doesn't coincide with the cells, we must set up the source terms before solving the PDE.
    // Pass in already updated CellPdeElementMap to speed up finding cells.
    this->SetUpSourceTermsForAveragedSourcePde(this->mpFeMesh, &this->mCellPdeElementMap);

    // Use SimpleLinearParabolicSolver as averaged Source PDE
    SimpleLinearParabolicSolver<DIM,DIM> solver(this->mpFeMesh,
                                                boost::static_pointer_cast<AbstractLinearParabolicPde<DIM,DIM> >(this->GetPde()).get(),
                                                p_bcc.get());

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
    AbstractBoxDomainPdeModifier<DIM>::SetupSolve(rCellPopulation,outputDirectory);

    // Copy the cell data to mSolution (this is the initial condition)
    SetupInitialSolutionVector(rCellPopulation);

    // Output the initial conditions on FeMesh
    this->UpdateAtEndOfOutputTimeStep(rCellPopulation);
}

template<unsigned DIM>
std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > MolecularConcentrationsDomainPdeModifier<DIM>::ConstructBoundaryConditionsContainer(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc(new BoundaryConditionsContainer<DIM,DIM,1>(false));

    if (!this->mSetBcsOnBoxBoundary)
    {
        EXCEPTION("Boundary conditions cannot yet be set on the cell population boundary for a MolecularConcentrationsDomainPdeModifier");
    }
    else // Apply BC at boundary nodes of box domain FE mesh
    {
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
    }

    return p_bcc;
}

template<unsigned DIM>
void MolecularConcentrationsDomainPdeModifier<DIM>::SetupInitialSolutionVector(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    double init_cond_value = mInitialValue;
    std::vector<double> init_cond(this->mpFeMesh->GetNumNodes());

    // TIP CELL INITIAL CONDITION TO MODIFY 
    // double TipCellXCoord = 2.0;
    // double TipCellYCoord = 0.0;
    // double TipCellZCoord = 0.0;

    // // initial coordinates of the tip cell (we suppose that there is only one tip cell at t=0 for the moment)
    // for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin(); cell_iter != rCellPopulation.End(); ++cell_iter)
    // {
    //     unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
    //     CellPtr pCell = rCellPopulation.GetCellUsingLocationIndex(node_index);

    //     if (pCell->GetCellProliferativeType()->IsType<TipCellMutationState>()){
    //         TipCellXCoord = rCellPopulation.GetLocationOfCellCentre(pCell)[0];
    //         TipCellYCoord = rCellPopulation.GetLocationOfCellCentre(pCell)[1];
    //         TipCellZCoord = rCellPopulation.GetLocationOfCellCentre(pCell)[2];
    //     }
    // }

    // PRINT_VARIABLE(TipCellXCoord)
    // PRINT_VARIABLE(TipCellYCoord)
    // PRINT_VARIABLE(TipCellZCoord)

    // initial coordinates of the endometriotic lesion 
    for(unsigned i=0; i<this->mpFeMesh->GetNumNodes(); i++){
        if(this->mpFeMesh->GetNode(i)->rGetLocation()[0] == mBoundaryCuboidMax){
            init_cond[i] = init_cond_value;
        } else {
            init_cond[i] = 0.0;
        }

        // TIP CELL INITIAL CONDITION TO MODIFY
        // if(this->mpFeMesh->GetNode(i)->rGetLocation()[0] == TipCellXCoord && this->mpFeMesh->GetNode(i)->rGetLocation()[1] == TipCellYCoord && this->mpFeMesh->GetNode(i)->rGetLocation()[2] == TipCellZCoord){
        //     init_cond[i] = 0.1*init_cond_value;
        // }
    }

    // Initialise mSolution
    this->mSolution = PetscTools::CreateVec(init_cond);

}

template<unsigned DIM>
void MolecularConcentrationsDomainPdeModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractBoxDomainPdeModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class MolecularConcentrationsDomainPdeModifier<1>;
template class MolecularConcentrationsDomainPdeModifier<2>;
template class MolecularConcentrationsDomainPdeModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MolecularConcentrationsDomainPdeModifier)