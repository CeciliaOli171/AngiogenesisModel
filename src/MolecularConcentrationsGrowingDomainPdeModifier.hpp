#ifndef MOLECULARCONCENTRATIONSGROWINGDOMAINPDEMODIFIER_HPP_
#define MOLECULARCONCENTRATIONSGROWINGDOMAINPDEMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractGrowingDomainPdeModifier.hpp"
#include "ParabolicGrowingDomainPdeModifier.hpp"
#include "BoundaryConditionsContainer.hpp"

/**
 * Subclass of ParabolicGrowingDomainPdeModifier.
 * A modifier class in which the PDE of VEGF concentration coupled to a cell-based simulation for angiogenesis is solved on a growing domain. 
 *
 * The value of the dependent variable at each cell is stored and updated in a CellData item.
 *
 * At each time step, the finite element mesh used to solve the PDE numerically
 * is defined by the spatial domain associated with the cell population. The
 * precise definition of this domain is implemented in the method
 * GetTetrahedralMeshForPdeModifier(), which is overridden for each cell population
 * class and is used in the AbstractGrowingDomainPdeModifier method GenerateFeMesh()
 * that is inherited by this class.
 */

template<unsigned DIM>
class MolecularConcentrationsGrowingDomainPdeModifier : public ParabolicGrowingDomainPdeModifier<DIM>
{
    friend class TestForcesModel;
    friend class TestAngiogenesisModelWithVegfConcentrationPde;

private:
    /* parameters */
    double mConstantBackground;

    /* serialisation */
    friend class boost::serialization::access;

    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractGrowingDomainPdeModifier<DIM> >(*this);
        archive & mConstantBackground;
    }

public:

    /**
     * Constructor.
     *
     * @param pPde A shared pointer to a linear PDE object (defaults to NULL)
     * @param pBoundaryCondition A shared pointer to an abstract boundary condition (defaults to NULL, corresponding to a constant boundary condition with value zero)
     * @param isNeumannBoundaryCondition Whether the boundary condition is Neumann (defaults to true)
     * @param solution solution vector (defaults to NULL)
     * @param initialValue initial value of vegf concentration
     * @param constantBackground constant baseline of vegf concentration
     */
    MolecularConcentrationsGrowingDomainPdeModifier(boost::shared_ptr<AbstractLinearPde<DIM,DIM> > pPde=boost::shared_ptr<AbstractLinearPde<DIM,DIM> >(), boost::shared_ptr<AbstractBoundaryCondition<DIM> > pBoundaryCondition=boost::shared_ptr<AbstractBoundaryCondition<DIM> >(), bool isNeumannBoundaryCondition=true, Vec solution=nullptr, double constantBackground=0.1);

    /**
     * Destructor.
     */
    ~MolecularConcentrationsGrowingDomainPdeModifier();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specifies what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specifies what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Helper method to construct the boundary conditions container for the PDE.
     *
     * @return the full boundary conditions container
     */ 
    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > ConstructBoundaryConditionsContainer();

    /**
     * Helper method to copy the CellData to the PDE solution.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateSolutionVectorVEGF(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Helper method to initialise the CellData for VEGF concentration.
     *
     * @param rCellPopulation reference to the cell population
     */
    void SetInitialCellDataVEGF(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    // void GenerateFeMeshAngiogenesis(AbstractCellPopulation<DIM,DIM>& rCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MolecularConcentrationsGrowingDomainPdeModifier)

namespace boost
{
namespace serialization
{
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const MolecularConcentrationsGrowingDomainPdeModifier<DIM> * t, const unsigned int file_version)
{
    if (t->GetSolution())
    {
        std::string archive_filename = ArchiveLocationInfo::GetArchiveDirectory() + "solution.vec";
        PetscTools::DumpPetscObject(t->GetSolution(), archive_filename);
    }
}

template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, MolecularConcentrationsGrowingDomainPdeModifier<DIM> * t, const unsigned int file_version)
{
    Vec solution = nullptr;

    std::string archive_filename = ArchiveLocationInfo::GetArchiveDirectory() + "solution.vec";
    FileFinder file_finder(archive_filename, RelativeTo::Absolute);

    if (file_finder.Exists())
    {
        PetscTools::ReadPetscObject(solution, archive_filename);
    }

    ::new(t)MolecularConcentrationsGrowingDomainPdeModifier<DIM>(boost::shared_ptr<AbstractLinearPde<DIM, DIM> >(),
                                                   boost::shared_ptr<AbstractBoundaryCondition<DIM> >(),
                                                   true,
                                                   solution);
}
}
} // namespace ...

#endif /*MOLECULARCONCENTRATIONSGROWINGDOMAINPDEMODIFIER_HPP_*/