#ifndef VEGFEQUATIONCELLPDE_HPP_
#define VEGFEQUATIONCELLPDE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellPopulation.hpp"
#include "AbstractLinearParabolicPde.hpp"
#include "CellwiseSourceParabolicPde.hpp"

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
 *
 * For a node of the finite element mesh with location x, the function rho(x)
 * equals one if there is a vessel tip associated with x, and zero otherwise. 
 */

template<unsigned DIM>
class VegfEquationCellPde : public CellwiseSourceParabolicPde<DIM>
{
    friend class TestForcesModel;

private:
    /* parameters */
    AbstractCellPopulation<DIM, DIM>& mrCellPopulation;
    double mDuDtCoefficient;
    double mDiffusionCoefficient;
    double mDecayCoefficient;
    double mCreationCoefficient;
    double mConsumptionCoefficient;

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
     * Serialize the PDE and its member variables.
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
    VegfEquationCellPde(AbstractCellPopulation<DIM, DIM>& rCellPopulation, double duDtCoefficient=1.0, double diffusionCoefficient=1.0, double decayCoefficient=1.0, double creationCoefficient=0.1, double consumptionCoefficient=0.01);

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
     * Never called.
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
     * Computes the creation, decay and consumption by vessel tips terms.
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
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VegfEquationCellPde)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a VegfEquationCellPde.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const VegfEquationCellPde<DIM>* t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM, DIM>* p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise a VegfEquationCellPde.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, VegfEquationCellPde<DIM>* t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM, DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)VegfEquationCellPde<DIM>(*p_cell_population);
}
}
} // namespace ...

#endif /*VEGFEQUATIONCELLPDE_HPP_*/