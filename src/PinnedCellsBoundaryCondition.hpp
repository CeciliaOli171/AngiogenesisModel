#ifndef PINNEDCELLSBOUNDARYCONDITION_HPP_
#define PINNEDCELLSBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"
#include "Debug.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * A boundary class to pin selected cells to certain positions.
 * Cells are pinned to whatever location they have when the constructor function is called
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class PinnedCellsBoundaryCondition : public AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM>
{
friend class TestForcesModel;

private:

    /**
     * Map to store the pinned nodes and their locations
     * - key: node index
     * - value: location vector node is pinned to
    */
    std::map<unsigned, c_vector<double,SPACE_DIM> > mPinnedNodes;

    /** Needed for serialization. */
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
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive & mPinnedNodes;
    }

public:

    /**
     * Constructor.
     *
     * @param pinned_node_indices a vector containing all the indices of the nodes that are to be pinned
     */
    PinnedCellsBoundaryCondition(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation, std::vector<unsigned> pinned_node_indices);

    /**
     * Overridden ImposeBoundaryCondition() method.
     *
     * Apply the cell population boundary conditions.
     *
     * @param rOldLocations the node locations before any boundary conditions are applied
     */
    void ImposeBoundaryCondition(const std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> >& rOldLocations);

    /**
     * Overridden VerifyBoundaryCondition() method.
     * Verify the boundary conditions have been applied.
     * This is called after ImposeBoundaryCondition() to ensure the condition is still satisfied.
     *
     * @return whether the boundary conditions are satisfied.
     */
    bool VerifyBoundaryCondition();

    /**
     * Function to get the IDs of the pinned nodes, required for load_construct_data
    */
    const std::vector<unsigned> GetPinnedNodeIDs() const;

    /**
     * Overridden OutputCellPopulationBoundaryConditionParameters() method.
     * Output cell population boundary condition parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(PinnedCellsBoundaryCondition)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a PinnedCellsBoundaryCondition.
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const PinnedCellsBoundaryCondition<ELEMENT_DIM, SPACE_DIM>* t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;

    // Archive the node indices (locations will be stored in the population)
    std::vector<unsigned> pinned_nodes = t->GetPinnedNodeIDs();

    // Flag for the end of the node index archive
    for (std::vector<unsigned>::const_iterator it = pinned_nodes.begin(); it != pinned_nodes.end(); it++)
    {
        ar << (*it);
    }
    ar << -1;
}

/**
 * De-serialize constructor parameters and initialize a PinnedCellsBoundaryCondition.
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, PinnedCellsBoundaryCondition<ELEMENT_DIM, SPACE_DIM>* t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population;
    ar >> p_cell_population;
    
    // Get all the node indices
    std::vector<unsigned> node_indices;
    unsigned max_nodes = p_cell_population->GetNumNodes();
    for (unsigned i = 0; i < max_nodes; i++)
    {
        int index;
        ar >> index;
        if (index == -1)
        {
            break;
        } else if (i == (max_nodes-1))
        {
            EXCEPTION("Error reading in the node indices in load_construct_data for PinnedCellsBoundaryCondition, didn't detect the end of the node list.");
        } else
        {
            node_indices.push_back(index);
        }
    }

    // Invoke inplace constructor to initialise instance
    ::new(t)PinnedCellsBoundaryCondition<ELEMENT_DIM, SPACE_DIM>(p_cell_population, node_indices);
}
}
} // namespace ...

#endif /*PINNEDCELLSBOUNDARYCONDITION_HPP_*/