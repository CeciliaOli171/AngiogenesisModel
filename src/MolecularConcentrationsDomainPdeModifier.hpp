#ifndef MOLECULARCONCENTRATIONSDOMAINPDEMODIFIER_HPP_
#define MOLECULARCONCENTRATIONSDOMAINPDEMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractBoxDomainPdeModifier.hpp"
#include "BoundaryConditionsContainer.hpp"


template<unsigned DIM>
class MolecularConcentrationsDomainPdeModifier : public AbstractBoxDomainPdeModifier<DIM>
{
private:
    double mBoundaryCuboidMax;
    double mInitialValue;

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractBoxDomainPdeModifier<DIM> >(*this);
    }

public:

    // constructor
    MolecularConcentrationsDomainPdeModifier(boost::shared_ptr<AbstractLinearPde<DIM,DIM> > pPde=boost::shared_ptr<AbstractLinearPde<DIM,DIM> >(),
                                  boost::shared_ptr<AbstractBoundaryCondition<DIM> > pBoundaryCondition=boost::shared_ptr<AbstractBoundaryCondition<DIM> >(),
                                  bool isNeumannBoundaryCondition=true,
                                  boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid=boost::shared_ptr<ChasteCuboid<DIM> >(),
                                  double stepSize=1.0,
                                  Vec solution=nullptr, 
                                  double boundaryCuboidMax=20.0,
                                  double initialValue=0.1);

    // destructor
    ~MolecularConcentrationsDomainPdeModifier();

    void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > ConstructBoundaryConditionsContainer(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    void SetupInitialSolutionVector(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

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
                                               boost::shared_ptr<ChasteCuboid<DIM> >(),
                                               1.0,
                                               solution);
}
}
} 

#endif /*MOLECULARCONCENTRATIONSDOMAINPDEMODIFIER_HPP_*/