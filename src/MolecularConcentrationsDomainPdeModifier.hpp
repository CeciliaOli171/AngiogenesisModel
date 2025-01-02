#ifndef MOLECULARCONCENTRATIONSDOMAINPDEMODIFIER_HPP_
#define MOLECULARCONCENTRATIONSDOMAINPDEMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractGrowingDomainPdeModifier.hpp"
#include "BoundaryConditionsContainer.hpp"


template<unsigned DIM>
class MolecularConcentrationsDomainPdeModifier : public AbstractGrowingDomainPdeModifier<DIM>
{
private:

    friend class boost::serialization::access;
   
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractGrowingDomainPdeModifier<DIM> >(*this);
    }

public:

    /**
     * Constructor.
     *
     * @param pPde A shared pointer to a linear PDE object (defaults to NULL)
     * @param pBoundaryCondition A shared pointer to an abstract boundary condition
     *     (defaults to NULL, corresponding to a constant boundary condition with value zero)
     * @param isNeumannBoundaryCondition Whether the boundary condition is Neumann (defaults to true)
     * @param solution solution vector (defaults to NULL)
     */
    MolecularConcentrationsDomainPdeModifier(boost::shared_ptr<AbstractLinearPde<DIM,DIM> > pPde=boost::shared_ptr<AbstractLinearPde<DIM,DIM> >(),
                                      boost::shared_ptr<AbstractBoundaryCondition<DIM> > pBoundaryCondition=boost::shared_ptr<AbstractBoundaryCondition<DIM> >(),
                                      bool isNeumannBoundaryCondition=true,
                                      Vec solution=nullptr);

    /**
     * Destructor.
     */
    virtual ~MolecularConcentrationsDomainPdeModifier();

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
    void UpdateSolutionVector(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MolecularConcentrationsDomainPdeModifier)

namespace boost
{
namespace serialization
{
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const MolecularConcentrationsDomainPdeModifier<DIM> * t, const unsigned int file_version)
{
    if (t->GetSolution())
    {
        std::string archive_filename = ArchiveLocationInfo::GetArchiveDirectory() + "solution.vec";
        PetscTools::DumpPetscObject(t->GetSolution(), archive_filename);
    }
}

template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, MolecularConcentrationsDomainPdeModifier<DIM> * t, const unsigned int file_version)
{
    Vec solution = nullptr;

    std::string archive_filename = ArchiveLocationInfo::GetArchiveDirectory() + "solution.vec";
    FileFinder file_finder(archive_filename, RelativeTo::Absolute);

    if (file_finder.Exists())
    {
        PetscTools::ReadPetscObject(solution, archive_filename);
    }

    ::new(t)MolecularConcentrationsDomainPdeModifier<DIM>(boost::shared_ptr<AbstractLinearPde<DIM, DIM> >(),
                                                   boost::shared_ptr<AbstractBoundaryCondition<DIM> >(),
                                                   true,
                                                   solution);
}
}
} // namespace ...

#endif /*MOLECULARCONCENTRATIONSDOMAINPDEMODIFIER_HPP_*/
