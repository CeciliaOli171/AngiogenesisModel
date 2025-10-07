#ifndef DAUGHTERCELLMODIFIER_HPP_
#define DAUGHTERCELLMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/* Header Files */
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "ArchiveOpener.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "CellPropertyCollection.hpp"
#include "CellPropertyRegistry.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "AbstractCentreBasedDivisionRule.hpp"
#include "AbstractCellBasedSimulationModifier.hpp"

/**
 * A modifier class (to implement setup, update and finalise methods), for use in angiogenesis simulations.
 */

template<unsigned DIM>
class DaughterCellModifier : public AbstractCellBasedSimulationModifier<DIM, DIM>
{
    friend class TestForcesModel;
    friend class TestAngiogenesisModel;
    friend class TestAngiogenesisModelWithVegfConcentrationPde;
    friend class TestAngiogenesisModelWithVegfConcentrationConstant;
    friend class TestAngiogenesisModelWithVegfConcentrationAnalyticalApproximationOfPde;

private:
    /* parameters */
    unsigned mHighestBranch;
    unsigned mHighestLoop;
    double mAnastomosisLength;
    double mThresholdLength;

    /**
     * @return the minimal distance between cells to have an anastomosis event.
     */
    double GetAnastomosisLength();

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
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM,DIM> >(*this);
        archive & mHighestBranch;
        archive & mHighestLoop;
        archive & mAnastomosisLength;
        archive & mThresholdLength;
    }

public:
    /**
     * Constructor.
     *
     * @param anastomosisLength maximal length between cell to have an anastomosis event.
     * @param thresholdLength minimal branch length
     */
    DaughterCellModifier(double anastomosisLength=1.5, double thresholdLength=2.0);

    /**
     * Destructor.
     */
    ~DaughterCellModifier();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specify what to do in the simulation at the end of each timestep.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Computes the indices of the cell's neighbours that are under the anastomosis threshold length.
     * @param rCellPopulation  The centre-based cell population
     * @param p_node_population the node population
     * @param pParentCell  The cell to divide
     *
     * @return the set of indices of cells that can go through anastomosis with the cell considered.
     */
    std::set<unsigned> GetAnastomosisNeighbours(AbstractCellPopulation<DIM, DIM>& rCellPopulation,NodeBasedCellPopulation<DIM>* p_node_population, CellPtr pParentCell);

    /**
     * Checks if the cell positioned next to the parent cell is a branching segment or not.
     * 
     * @param rCellPopulation  The centre-based cell population
     * @param p_node_population the node population
     * @param pParentCell  The cell to divide
     *
     * @return a boolean depending on the cell's neighbour type.
     */
    bool IsBranchingSegmentNextToCell(AbstractCellPopulation<DIM, DIM>& rCellPopulation,NodeBasedCellPopulation<DIM>* p_node_population, CellPtr pParentCell);

    /**
     * Computes the vector of the position of the closest neighbour to the cell.
     *
     * @param rCellPopulation reference to the cell population
     * @param pCell the cell
     * @param neighbouring_node_indices a set of the indices of the neighbours of a node in the same branch
     *
     * @return the coordinates of the closest neighbour to the cell.
     */
    std::pair<c_vector<double, DIM>, unsigned> ClosestNeighbour(AbstractCellPopulation<DIM, DIM>& rCellPopulation, CellPtr pCell, std::set<unsigned> neighbouring_node_indices);

    /**
     * Computes the anastomosis event: updates Cell Data.
     *
     * @param rCellPopulation reference to the cell population
     * @param pCell the cell
     * @param neighbouring_node_indices a set of the indices of the neighbours of a node in the same branch
     */
    void CalculateAnastomosisVector(AbstractCellPopulation<DIM, DIM>& rCellPopulation, NodeBasedCellPopulation<DIM>* p_node_population, CellPtr pParentCell);

    /**
     * Overridden SetupSolve() method.
     *
     * Specify what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Updates cell data after anastomosis and cell division.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     *
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DaughterCellModifier)

#endif /*DAUGHTERCELLMODIFIER_HPP_*/
