#include "MolecularConcentrationsBoxDomainPdeModifier.hpp"
#include "SimpleLinearParabolicSolver.hpp"
#include "VesselTipMutationState.hpp"
#include "AbstractPdeModifier.hpp"

#include "Debug.hpp"

template<unsigned DIM>
MolecularConcentrationsBoxDomainPdeModifier<DIM>::MolecularConcentrationsBoxDomainPdeModifier(boost::shared_ptr<AbstractLinearPde<DIM,DIM> > pPde, boost::shared_ptr<AbstractBoundaryCondition<DIM> > pBoundaryCondition, bool isNeumannBoundaryCondition, boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid, double stepSize, Vec solution, double sourceLocation, double sourceValue, double constantBackground) : ParabolicBoxDomainPdeModifier<DIM>(pPde, pBoundaryCondition, isNeumannBoundaryCondition, pMeshCuboid, stepSize, solution), mSourceLocation(sourceLocation), mSourceValue(sourceValue), mConstantBackground(constantBackground)
{
}

template<unsigned DIM>
MolecularConcentrationsBoxDomainPdeModifier<DIM>::~MolecularConcentrationsBoxDomainPdeModifier()
{
}

template<unsigned DIM>
void MolecularConcentrationsBoxDomainPdeModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    AbstractBoxDomainPdeModifier<DIM>::SetupSolve(rCellPopulation,outputDirectory);

    // Copy the cell data to mSolution (this is the initial condition)
    SetupInitialSolutionVectorVEGF(rCellPopulation);

    // Output the initial conditions on FeMesh
    this->UpdateAtEndOfOutputTimeStep(rCellPopulation);
}

template<unsigned DIM>
std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > MolecularConcentrationsBoxDomainPdeModifier<DIM>::ConstructBoundaryConditionsContainer(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc(new BoundaryConditionsContainer<DIM,DIM,1>(false));

    if (!this->mSetBcsOnBoxBoundary)
    {
        EXCEPTION("Boundary conditions cannot yet be set on the cell population boundary for a MolecularConcentrationsDomainPdeModifier");
    }
    else 
    {
        // those loops could be merged (to do later for optimisation)
        // no-flux boundary condition
        for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator elem_iter = this->mpFeMesh->GetBoundaryElementIteratorBegin(); elem_iter != this->mpFeMesh->GetBoundaryElementIteratorEnd(); ++elem_iter)
        {
            if(((*elem_iter)->GetNodeLocation(0))[0] != mSourceLocation){
                p_bcc->AddNeumannBoundaryCondition(*elem_iter, this->mpBoundaryCondition.get());
            }
        }
        // source term at the source location inside the lesion
        for (typename TetrahedralMesh<DIM,DIM>::BoundaryNodeIterator node_iter = this->mpFeMesh->GetBoundaryNodeIteratorBegin(); node_iter != this->mpFeMesh->GetBoundaryNodeIteratorEnd(); ++node_iter)
        {
            if((*node_iter)->GetPoint()[0] == mSourceLocation){
                p_bcc->AddDirichletBoundaryCondition(*node_iter, this->mpBoundaryCondition.get());
            }
        }
    }

    return p_bcc;
}

template<unsigned DIM>
void MolecularConcentrationsBoxDomainPdeModifier<DIM>::SetupInitialSolutionVectorVEGF(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    double source_value = mSourceValue;
    double constant_background = mConstantBackground;
    std::vector<double> init_cond(this->mpFeMesh->GetNumNodes());

    // initial coordinates of the endometriotic lesion 
    for(unsigned i=0; i<this->mpFeMesh->GetNumNodes(); i++){
        if(this->mpFeMesh->GetNode(i)->rGetLocation()[0] == mSourceLocation){
            init_cond[i] = source_value;
        } else {
            init_cond[i] = constant_background;
        }
    }

    // Initialise mSolution
    this->mSolution = PetscTools::CreateVec(init_cond);
}

// Explicit instantiation
template class MolecularConcentrationsBoxDomainPdeModifier<1>;
template class MolecularConcentrationsBoxDomainPdeModifier<2>;
template class MolecularConcentrationsBoxDomainPdeModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MolecularConcentrationsBoxDomainPdeModifier)