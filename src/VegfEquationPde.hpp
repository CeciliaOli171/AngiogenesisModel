#ifndef VEGFEQUATIONPDE_HPP_
#define VEGFEQUATIONPDE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractLinearParabolicPde.hpp"
#include "AveragedSourceParabolicPde.hpp"

template<unsigned DIM>
class VegfEquationPde : public AveragedSourceParabolicPde<DIM>
{
    friend class TestCellBasedParabolicPdes;

private:

    friend class boost::serialization::access;

    AbstractCellPopulation<DIM, DIM>& mrCellPopulation;
    double mDuDtCoefficient;
    double mDiffusionCoefficient;
    double mDecayCoefficient;
    double mCreationCoefficient;
    double mConsumptionCoefficient;
    double mConstantBackground;
    std::vector<double> mCellDensityOnCoarseElements;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
       archive & boost::serialization::base_object<AbstractLinearParabolicPde<DIM, DIM> >(*this);
       archive & mDuDtCoefficient;
       archive & mDiffusionCoefficient;
       archive & mDecayCoefficient;
       archive & mCreationCoefficient;
       archive & mConsumptionCoefficient;
       archive & mConstantBackground;
       archive & mCellDensityOnCoarseElements;
    }

public:

    VegfEquationPde(AbstractCellPopulation<DIM, DIM>& rCellPopulation, double duDtCoefficient=1.0, double diffusionCoefficient=1.0, double decayCoefficient=1.0, double creationCoefficient=0.1, double consumptionCoefficient=0.01, double constantBackground=0.1);

    void SetupSourceTerms(TetrahedralMesh<DIM,DIM>& rCoarseMesh, std::map<CellPtr, unsigned>* pCellPdeElementMap=nullptr);

    double ComputeDuDtCoefficientFunction(const ChastePoint<DIM>& rX);

    double ComputeSourceTerm(const ChastePoint<DIM>& rX, double u, Element<DIM,DIM>* pElement=NULL);

    double ComputeSourceTermAtNode(const Node<DIM>& rNode, double u);

    c_matrix<double,DIM,DIM> ComputeDiffusionTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement=NULL);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VegfEquationPde)

namespace boost
{
namespace serialization
{
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const VegfEquationPde<DIM>* t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM, DIM>* p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, VegfEquationPde<DIM>* t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM, DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)AveragedSourceParabolicPde<DIM>(*p_cell_population);
}
}
} 

#endif /*VEGFEQUATIONPDE_HPP_*/
