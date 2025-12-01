#ifndef MOLECULARCONCENTRATIONSBOXDOMAINPDEMODIFIER_HPP_
#define MOLECULARCONCENTRATIONSBOXDOMAINPDEMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractBoxDomainPdeModifier.hpp"
#include "ParabolicBoxDomainPdeModifier.hpp"
#include "BoundaryConditionsContainer.hpp"

/**
 * Subclass of ParabolicBoxDomainPdeModifier.
 * A modifier class in which the PDE of VEGF concentration coupled to a cell-based simulation for angiogenesis is solved on a coarse domain.
 *
 * The finite element mesh used to solve the PDE numerically is a fixed tessellation of
 * a cuboid (box), which must be supplied to the constructor. The value of the dependent
 * variable is interpolated between coarse mesh nodes to obtain a value at each cell,
 * which is stored and updated in a CellData item.
 *
 * At each time step the boundary condition supplied to the constructor may be imposed
 * either on the boundary of the box domain, or on the boundary of the cell population
 * (which is assumed to lie within the box domain). This choice can be made using the
 * AbstractBoxDomainPdeModifier method SetBcsOnBoxBoundary(), which is inherited by this
 * class.
 */


template<unsigned DIM>
class MolecularConcentrationsBoxDomainPdeModifier : public ParabolicBoxDomainPdeModifier<DIM>
{
    friend class TestForcesModel;
    friend class TestAngiogenesisModelWithVegfConcentrationPde;

private:
    /* parameters */
    double mSourceLocation;
    double mSourceValue;
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
        archive & boost::serialization::base_object<AbstractBoxDomainPdeModifier<DIM> >(*this);
        archive & mSourceLocation;
        archive & mSourceValue;
        archive & mConstantBackground;
    }

public:
    /**
     * Constructor.
     *
     * @param pPde A shared pointer to a linear PDE object (defaults to NULL)
     * @param pBoundaryCondition A shared pointer to an abstract boundary condition (defaults to NULL, corresponding to a constant boundary condition with value zero)
     * @param isNeumannBoundaryCondition Whether the boundary condition is Neumann (defaults to true)
     * @param pMeshCuboid A shared pointer to a ChasteCuboid specifying the outer boundary for the FE mesh (defaults to NULL)
     * @param stepSize step size to be used in the FE mesh (defaults to 1.0, i.e. the default cell size)
     * @param solution solution vector (defaults to NULL)
     * @param sourceLocation location of source inside lesion
     * @param sourceValue value of vegf concentration at lesion source
     * @param constantBackground constant baseline of vegf concentration
     */
    MolecularConcentrationsBoxDomainPdeModifier(boost::shared_ptr<AbstractLinearPde<DIM,DIM> > pPde=boost::shared_ptr<AbstractLinearPde<DIM,DIM> >(),boost::shared_ptr<AbstractBoundaryCondition<DIM> > pBoundaryCondition=boost::shared_ptr<AbstractBoundaryCondition<DIM> >(), bool isNeumannBoundaryCondition=true, boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid=boost::shared_ptr<ChasteCuboid<DIM> >(), double stepSize=1.0, Vec solution=nullptr, double sourceLocation=0.0, double sourceValue=1.0, double constantBackground=0.1);

    /**
     * Destructor.
     */
    ~MolecularConcentrationsBoxDomainPdeModifier();

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
     * Overridden ConstructBoundaryConditionsContainer() method.
     *
     * Helper method to construct the boundary conditions container for the PDE.
     *
     * @param rCellPopulation reference to the cell population
     *
     * @return the full boundary conditions container
     */
    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > ConstructBoundaryConditionsContainer(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Helper method to initialise the PDE solution using the CellData.
     *
     * Here we assume for the initial condition that the 
     *
     * @param rCellPopulation reference to the cell population
     */
    void SetupInitialSolutionVectorVEGF(AbstractCellPopulation<DIM,DIM>& rCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MolecularConcentrationsBoxDomainPdeModifier)

namespace boost
{
namespace serialization
{
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const MolecularConcentrationsBoxDomainPdeModifier<DIM> * t, const unsigned int file_version)
{
    if (t->GetSolution())
    {
        std::string archive_filename = ArchiveLocationInfo::GetArchiveDirectory() + "solution.vec";
        PetscTools::DumpPetscObject(t->GetSolution(), archive_filename);
    }
}

template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, MolecularConcentrationsBoxDomainPdeModifier<DIM> * t, const unsigned int file_version)
{
    Vec solution = nullptr;

    std::string archive_filename = ArchiveLocationInfo::GetArchiveDirectory() + "solution.vec";
    FileFinder file_finder(archive_filename, RelativeTo::Absolute);

    if (file_finder.Exists())
    {
        PetscTools::ReadPetscObject(solution, archive_filename);
    }

    ::new(t)MolecularConcentrationsBoxDomainPdeModifier<DIM>(boost::shared_ptr<AbstractLinearPde<DIM, DIM> >(),boost::shared_ptr<AbstractBoundaryCondition<DIM> >(), true, boost::shared_ptr<ChasteCuboid<DIM> >(), 1.0, solution, 0.0, 1.0, 0.1);
}
}
} 

#endif /*MOLECULARCONCENTRATIONSDOMAINPDEMODIFIER_HPP_*/