#ifndef VEGFEQUATIONPDE_HPP_
#define VEGFEQUATIONPDE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractLinearParabolicPde.hpp"
#include "AveragedSourceParabolicPde.hpp"

/**
 * A parabolic reaction-diffusion PDE of vegf concentration to be solved numerically using the finite element method, 
 * for coupling to a cell-based simulation for angiogenesis.
 *
 * The PDE takes the form:
 * c*du/dt = Grad.(Dc*Grad(u)) + (Ac-Mc)u - Epsc*n(x)*u,
 *
 * where the scalars c, Dc, Ac, Mc and Epsc are specified by the members mDuDtCoefficient,
 * mDiffusionCoefficient, mCreationCoefficient, mDecayCoefficient and mConsumptionCoefficient, respectively. 
 * Their values must be set in the constructor. The function n(x) denotes the local density of non-apoptotic cells. 
 * This quantity is computed for each element of a 'coarse' finite element mesh that is
 * passed to the method SetupSourceTerms() and stored in the member mCellDensityOnCoarseElements.
 * For a point x, rho(x) is defined to be the number of non-apoptotic cells whose
 * centres lie in each finite element containing that point, scaled by the area of
 * that element.
 */

template<unsigned DIM>
class VegfEquationPde : public AveragedSourceParabolicPde<DIM>
{
    friend class TestCellBasedParabolicPdes;
    friend class TestForcesModel;

private:
    /* parameters */
    AbstractCellPopulation<DIM, DIM>& mrCellPopulation;
    double mDuDtCoefficient;
    double mDiffusionCoefficient;
    double mDecayCoefficient;
    double mCreationCoefficient;
    double mConsumptionCoefficient;
    std::vector<double> mCellDensityOnCoarseElements;

    /**
     * @return vegf rate of change term coefficient.
     */
    double GetDuDtCoefficient();

    /**
     * @return vegf diffusion coefficient.
     */
    double GetDiffusionCoefficient();

    /**
     * @return vegf decay coefficient.
     */
    double GetDecayCoefficient();

    /**
     * @return vegf creation coefficient.
     */
    double GetCreationCoefficient();

    /**
     * @return vegf consumption by vessel tips coefficient.
     */
    double GetConsumptionCoefficient();

    /* serialisation */
    friend class boost::serialization::access;

    /**
     * Serialize the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
       archive & boost::serialization::base_object<AbstractLinearParabolicPde<DIM, DIM> >(*this);
       archive & mDuDtCoefficient;
       archive & mDiffusionCoefficient;
       archive & mDecayCoefficient;
       archive & mCreationCoefficient;
       archive & mConsumptionCoefficient;
       archive & mCellDensityOnCoarseElements;
    }

public:
    /**
     * Constructor.
     *
     * @param rCellPopulation reference to the cell population
     * @param duDtCoefficient rate of reaction 
     * @param diffusionCoefficient rate of diffusion 
     * @param decayCoefficient rate of decay 
     * @param creactionCoefficient rate of creation 
     * @param consumptionCoefficient rate of consumption by vessel tips
     */
    VegfEquationPde(AbstractCellPopulation<DIM, DIM>& rCellPopulation, double duDtCoefficient=1.0, double diffusionCoefficient=1.0, double decayCoefficient=1.0, double creationCoefficient=0.1, double consumptionCoefficient=0.01);

    /**
     * @return const reference to the cell population (used in archiving).
     */
    const AbstractCellPopulation<DIM>& rGetCellPopulation() const;

    /**
     * Overridden SetupSourceTerms() method. 
     *
     * Set up the source terms.
     *
     * @param rCoarseMesh reference to the coarse mesh
     * @param pCellPdeElementMap optional pointer to the map from cells to coarse elements
     */
    void SetupSourceTerms(TetrahedralMesh<DIM,DIM>& rCoarseMesh, std::map<CellPtr, unsigned>* pCellPdeElementMap=nullptr);

    /**
     * Overridden ComputeDuDtCoefficientFunction() method.
     *
     * Computes the reaction rate term.
     *
     * @param rX the point in space at which the function c is computed
     *
     * @return the function c(x) in "c(x) du/dt = Grad.(DiffusionTerm(x)*Grad(u))+LinearSourceTerm(x)+NonlinearSourceTerm(x, u)"
     */
    double ComputeDuDtCoefficientFunction(const ChastePoint<DIM>& rX);

    /**
     * Overridden ComputeSourceTerm() method.
     *
     * Computes the creation, decay and consumption by vessel tips terms.
     *
     * @param rX the point in space at which the nonlinear source term is computed
     * @param u the value of the dependent variable at the point
     * @param pElement the mesh element that x is contained in (optional; defaults to NULL).
     *
     * @return computed source term.
     */
    double ComputeSourceTerm(const ChastePoint<DIM>& rX, double u, Element<DIM,DIM>* pElement=NULL);

    /**
     * Overridden ComputeSourceTermAtNode() method. 
     *
     * Never called.
     *
     * @param rNode the node at which the nonlinear source term is computed
     * @param u the value of the dependent variable at the node
     *
     * @return computed source term at a node.
     */
    double ComputeSourceTermAtNode(const Node<DIM>& rNode, double u);

    /**
     * Overridden ComputeDiffusionTerm() method.
     *
     * Computes the diffusion term.
     *
     * @param rX the point in space at which the diffusion term is computed
     * @param pElement the mesh element that x is contained in (optional; defaults to NULL).
     *
     * @return a matrix.
     */
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
