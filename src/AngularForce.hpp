#ifndef ANGULARFORCE_HPP_
#define ANGULARFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "VesselSegmentMutationState.hpp"
#include "VesselTipMutationState.hpp"

/**
 * An angular force class. Describes the interactions between the vessel segment and its micro-environment that induce the alignment of the vessel segments in a vascular branch.
 */

template<unsigned DIM>
class AngularForce  : public AbstractForce<DIM>
{
friend class TestForcesModel;
friend class TestAngiogenesisModel;
friend class TestAngiogenesisModelWithVegfConcentrationPde;
friend class TestAngiogenesisModelWithVegfConcentrationConstant;
friend class TestAngiogenesisModelWithVegfConcentrationAnalyticalApproximationOfPde;

private:
    /* parameters */
    double mOmegaa;

    /**
     * @return the angular persistence spring constant.
     */
    double GetAngularPersistence();  

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
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mOmegaa;
    }

public:
    /**
     * Constructor.
     *
     * @param omegaa the angular persistence spring constant
     */
    AngularForce(double omegaa = 1E-5); 

    /**
     * Destructor.
     */ 
    ~AngularForce();

    /** 
     * Calculates the nth value of a neighbouring node indices set.
     * 
     * @param neighbouring_node_indices a set of the indices of the neighbours of a node in the same branch
     * @param n an integer
     *
     * @return the nth indice of the node's neighbours set
     */
    unsigned GetNthNeighbourIndice(std::set<unsigned> neighbouring_node_indices, int n);

    /**
     * Selects all the cells in the same branch.
     *
     * @param rCellPopulation reference to the cell population
     * @param p_node_population reference to the node population 
     * @param pCell the cell
     * 
     * @return a set with the indices of the neighbours of the cell that are in the same branch 
     */
    std::set<unsigned> GetSameBranchNeighbours(AbstractCellPopulation<DIM, DIM>& rCellPopulation,NodeBasedCellPopulation<DIM>* p_node_population, CellPtr pCell);

    /**
     * Calculates the angle between two vectors. 
     *
     * @param u a first vector 
     * @param v a second vector
     *
     * @return the angle between vectors u and v 
     */
    double GetAngleFromVectors(c_vector<double,DIM> u, c_vector<double,DIM> v);

    /**
     * Calculates the angle made between a cell and its two closest neighbours (branching point not considered).
     *
     * @param rCellPopulation reference to the cell population
     * @param pCell the cell
     * @param neighbouring_node_indices a set of the indices of the neighbours of a node in the same branch
     *
     * @return the angle and vector coordinates of the two closest neighbours
     */
    std::tuple<double, c_vector<double,DIM>, c_vector<double,DIM>> ClosestAngleVesselSegment(AbstractCellPopulation<DIM>& rCellPopulation, CellPtr pCell, std::set<unsigned> neighbouring_node_indices);

    /**
     * Calculates the smallest angle made between a cell and two of its neighbours (branching point not considered).
     *
     * @param rCellPopulation reference to the cell population
     * @param pCell the cell
     * @param neighbouring_node_indices a set of the indices of the neighbours of a node in the same branch
     *
     * @return the angle and vector coordinates of the two neighbours
     */
    std::tuple<double, c_vector<double,DIM>, c_vector<double,DIM>> OptimalAngleVesselElement(AbstractCellPopulation<DIM>& rCellPopulation, CellPtr pCell, std::set<unsigned> neighbouring_node_indices); 

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the cell population
     *
     * Fa = omegaa (alphaangular - pi) ((x_b - x_s) + (x_c - x_s))/(|(x_b - x_s) + (x_c - x_s)|)
     * (s a vessel segment, b and c its neighbours that optimise the value of alphaangular)
     * inspired by Perfhal et al. (2017)
     *
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AngularForce) 

#endif /*ANGULARFORCE_HPP_*/