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

template<unsigned DIM>
class DaughterCellModifier : public AbstractCellBasedSimulationModifier<DIM, DIM>
{
private:

    unsigned mHighestBranch;
    unsigned mHighestLoop;

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM,DIM> >(*this);
    }

public:

    DaughterCellModifier();

    ~DaughterCellModifier();

    void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    // calculates the neighbours in the anastomosis cut-off length that belong to a different branch
    std::set<unsigned> GetAnastomosisNeighbours(AbstractCellPopulation<DIM, DIM>& rCellPopulation,NodeBasedCellPopulation<DIM>* p_node_population, CellPtr pParentCell);

    // calculates the coordinates of the closest neighbour of a cell 
    std::pair<c_vector<double, DIM>, unsigned> ClosestNeighbour(AbstractCellPopulation<DIM, DIM>& rCellPopulation, CellPtr pCell, std::set<unsigned> neighbouring_node_indices);

    // calculates the daughter and parent position for a tip cell division in the case of anastomosis
    void CalculateAnastomosisVector(AbstractCellPopulation<DIM, DIM>& rCellPopulation, NodeBasedCellPopulation<DIM>* p_node_population, CellPtr pParentCell);

    void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DaughterCellModifier)

#endif /*DAUGHTERCELLMODIFIER_HPP_*/
