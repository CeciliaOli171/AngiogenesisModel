#include "PinnedCellsBoundaryCondition.hpp"
#include "NodeBasedCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PinnedCellsBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::PinnedCellsBoundaryCondition(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>* pCellPopulation,
                                                    std::vector<unsigned> pinned_node_indices)
        : AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM>(pCellPopulation)
{
    // Check that the pinned node indices isn't empty
    assert(!pinned_node_indices.empty());

    NodeBasedCellPopulation<SPACE_DIM>* pNodeBasedPopulation = dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(pCellPopulation);
    if (!pNodeBasedPopulation)
    {
        EXCEPTION("PinnedCellsBoundaryCondition only implemented and tested for Node Based Cell Populations.");
    }

    // Store the cell's initial (and hence pinned) locations
    for (std::vector<unsigned>::const_iterator it = pinned_node_indices.begin(); it != pinned_node_indices.end(); ++it)
    {
        // Get the node's location
        Node<SPACE_DIM>* pNode = pNodeBasedPopulation->GetNode(*it);
        c_vector<double,SPACE_DIM> node_location = pNode->rGetLocation();
        // Add it to the map
        mPinnedNodes[*it]=node_location;
    }
}  

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PinnedCellsBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::ImposeBoundaryCondition(const std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> >& rOldLocations)
{
    for (typename std::map<unsigned, c_vector<double,SPACE_DIM> >::const_iterator node_it = mPinnedNodes.begin(); 
          node_it != mPinnedNodes.end(); ++node_it)
    {
        // Get the node index and location from the stored values
        unsigned index = node_it->first;
        c_vector<double,SPACE_DIM> location = node_it->second;
        // Update the node's location to its pinned location
        Node<SPACE_DIM>* p_node = this->mpCellPopulation->GetNode(index);
        p_node->rGetModifiableLocation() = location;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool PinnedCellsBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::VerifyBoundaryCondition()
{
    // Loop over the stored pinned nodes and locations and check against the actual node locations
    for (typename std::map<unsigned, c_vector<double,SPACE_DIM> >::const_iterator node_it = mPinnedNodes.begin(); 
          node_it != mPinnedNodes.end(); ++node_it)
    {
        // Get the node index and location from the stored values
        unsigned index = node_it->first;
        c_vector<double,SPACE_DIM> stored_location = node_it->second;
        // Get the current location
        c_vector<double,SPACE_DIM> actual_location = this->mpCellPopulation->GetNode(index)->rGetLocation();
        // Check the node is in the correct location
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            if (actual_location[i] != stored_location[i])
            {
                return(false);
            }
        }
    }

    return(true);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const std::vector<unsigned> PinnedCellsBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::GetPinnedNodeIDs() const
{
   std::vector<unsigned> node_ids;
    for (typename std::map<unsigned, c_vector<double,SPACE_DIM> >::const_iterator node_it = mPinnedNodes.begin(); 
          node_it != mPinnedNodes.end(); ++node_it)
    {
        unsigned index = node_it->first;
        node_ids.push_back(index);
    }
   return(node_ids);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PinnedCellsBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<PinnedCells>\n";
    // Loop over the stored pinned nodes and locations and check against the actual node locations
    for (typename std::map<unsigned, c_vector<double,SPACE_DIM> >::const_iterator node_it = mPinnedNodes.begin(); 
          node_it != mPinnedNodes.end(); ++node_it)
    {
        // Get the node index and location from the stored values
        unsigned index = node_it->first;
        c_vector<double,SPACE_DIM> location = node_it->second;
        *rParamsFile << "\t\t\t<Cell>";
        *rParamsFile << index << ",";
        for (unsigned index=0; index != SPACE_DIM-1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
        {
            *rParamsFile << location[index] << ",";
        }
        *rParamsFile << location[SPACE_DIM-1] << "</Cell>\n";
    }
    *rParamsFile << "</PinnedCells>\n";

    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Explicit instantiation
template class PinnedCellsBoundaryCondition<1,1>;
template class PinnedCellsBoundaryCondition<1,2>;
template class PinnedCellsBoundaryCondition<2,2>;
template class PinnedCellsBoundaryCondition<1,3>;
template class PinnedCellsBoundaryCondition<2,3>;
template class PinnedCellsBoundaryCondition<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(PinnedCellsBoundaryCondition)