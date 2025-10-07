#ifndef SPROUTINGRULE_HPP_
#define SPROUTINGRULE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/* Header Files */
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "ArchiveOpener.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "CellPropertyCollection.hpp"
#include "CellPropertyRegistry.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "AbstractCentreBasedDivisionRule.hpp"
#include "ClassIsAbstract.hpp"
#include "CellData.hpp"

/**
 * A cell division rule class for use in centre-based simulations implemented for angiogenesis. The purpose of this class is to return the locations of the two daughters of a dividing cell depending on the type of division. Vessel tips can either divide into either a VS and a VT, lengthening the vessel (asymmetric division); or into two VTs, creating a new sprout (symmetric division).
 */

// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCentreBasedCellPopulation;
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCentreBasedDivisionRule;

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class SproutingRule  : public AbstractCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>
{
    friend class TestForcesModel;
    friend class TestAngiogenesisModel;
    friend class TestAngiogenesisModelWithVegfConcentrationPde;
    friend class TestAngiogenesisModelWithVegfConcentrationConstant;
    friend class TestAngiogenesisModelWithVegfConcentrationAnalyticalApproximationOfPde;

private:
    /* parameters */
    double mMaxSproutingRate;
    double mThresholdLength;

    /**
     * @return maximum sprouting rate
     */
    double GetMaxSproutingRate();

    /**
     * @return minimal branch length
     */
    double GetThresholdLength();

    /* serialisation */
    friend class boost::serialization::access;

    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // not in AbstractCentreBasedDivisionRule.hpp
        archive & boost::serialization::base_object<AbstractCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive & mMaxSproutingRate;
        archive & mThresholdLength;
    }

public:
    /**
     * Constructor.
     * 
     * @param maxSproutingRate the maximum sprouting rate 
     * @param thresholdLength the minimum branch length
     */
    SproutingRule(double maxSproutingRate = 0.08, double thresholdLength=2.0);

    /**
     * Destructor.
     */
    ~SproutingRule();

    /**
     * Computes the probability for a vessel tip to create a new sprout.
     * 
     * @param rCellPopulation the centre-based cell population
     * @param pParentCell the cell to divide
     *
     * @return sprouting probability.
     */
    virtual double GetSproutingProbability(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pParentCell);

    /**
     * Computes the set of all the vessel segment neighbours of the vessel tip.
     *
     * @param rCellPopulation  The centre-based cell population
     * @param p_node_population the node population
     * @param pParentCell  The cell to divide
     *
     * @return a set with the indices of the neighbours of the cell that are in the same branch and vessel segment
     */
    std::set<unsigned> GetVesselSegmentNeighbours(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, NodeBasedCellPopulation<SPACE_DIM>* p_node_population, CellPtr pParentCell); 

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
     * Computes the vector of the position of the closest neighbour to the cell.
     *
     * @param rCellPopulation reference to the cell population
     * @param pCell the cell
     * @param neighbouring_node_indices a set of the indices of the neighbours of a node in the same branch
     *
     * @return the coordinates of the closest neighbour to the cell.
     */
    std::pair<c_vector<double, SPACE_DIM>, unsigned> ClosestNeighbour(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pCell, std::set<unsigned> neighbouring_node_indices);

    /**
     * Checks if the cell positioned next to the parent cell is a branching segment or not.
     * 
     * @param rCellPopulation  The centre-based cell population
     * @param p_node_population the node population
     * @param pParentCell  The cell to divide
     *
     * @return a boolean depending on the cell's neighbour type.
     */
    bool IsBranchingSegmentNextToCell(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation,NodeBasedCellPopulation<SPACE_DIM>* p_node_population, CellPtr pParentCell);

    /**
     * Computes the position of the vector perpendicular to the vector of daughter cell.
     *
     * @param daughter_direction vector of daughter cell 
     * @param closest_neighbour vector of closest neighbour to the cell
     *
     * @return the coordinates of the perpendicular vector to daughter cell.
     */
    c_vector<double, SPACE_DIM> PerpendicularDaughterDirection(c_vector<double, SPACE_DIM> daughter_direction, c_vector<double, SPACE_DIM> closest_neighbour);

    /**
     * Return a pair of vectors used to position the daughters of a dividing vessel segment.
     * (not used for the moment)
     *
     * @param pParentCell  The cell to divide
     * @param rCellPopulation  The centre-based cell population
     *
     * @return the two daughter cell positions.
     */
    std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > CalculateCellDivisionVesselSegmentVector(CellPtr pParentCell, AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);

    /**
     * Computes a pair of vectors used to position the daughter cells obtained after a vessel tip division depending on the sprouting probability.
     *
     * @param rCellPopulation  The centre-based cell population
     * @param p_node_population the node population
     * @param pParentCell  The cell to divide
     *
     * @return the two daughter cell positions.
     */
    std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > CalculateGrowthVector( AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, NodeBasedCellPopulation<SPACE_DIM>* p_node_population, CellPtr pParentCell);

    /**
     * Return a pair of vectors used to position the daughters of a dividing vessel tip.
     *
     * @param pParentCell  The cell to divide
     * @param rCellPopulation  The centre-based cell population
     *
     * @return the two daughter cell positions.
     */
    std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > CalculateCellDivisionVesselTipVector(CellPtr pParentCell, AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);

    /**
     * Overridden CalculateCellDivisionVector() method.
     *
     * Return a pair of vectors used to position the daughters of a dividing cell.
     *
     * @param pParentCell  The cell to divide
     * @param rCellPopulation  The centre-based cell population
     *
     * @return the two daughter cell positions.
     */
    std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > CalculateCellDivisionVector(CellPtr pParentCell, AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);

    /**
     * Output the name of the concrete class and call OutputCellCentreBasedDivisionRuleParameters().
     *
     * @param rParamsFile  The stream of the parameter file
     */
    void OutputCellCentreBasedDivisionRuleInfo(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SproutingRule)

#endif /*SPROUTINGRULE_HPP_*/