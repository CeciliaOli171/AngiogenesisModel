#ifndef SPROUTING_HPP_
#define SPROUTING_HPP_

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
#include "AbstractCentreBasedDivisionRule.hpp"
#include "ClassIsAbstract.hpp"

// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCentreBasedCellPopulation;
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCentreBasedDivisionRule;

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class SproutingRule  : public AbstractCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>
{

private:

    c_vector<double, SPACE_DIM> mDaughterLocation; // we do not use it : remove it ? 

    // allow to archive the force model object in a cell-based simulation 
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // not in AbstractCentreBasedDivisionRule.hpp
        archive & boost::serialization::base_object<AbstractCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    // constructor 
    SproutingRule();

    // destructor 
    ~SproutingRule();

    // length of the vessel 
    double LengthVesselSegment(std::set<unsigned> neighbouring_node_indices, CellPtr cell_ptr, AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);

    // optimal angle between a vessel segment and its neighbours 
    double OptimalAngleVesselSegment(std::set<unsigned> neighbouring_node_indices, CellPtr cell_ptr, AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);

    // sprouting probability for the daughter cell 
    void DaughterTypeOfCell(AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rOldCellPopulation);

    // overrides CalculateCellDivisionVector
    std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > CalculateCellDivisionVector(CellPtr pParentCell, AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);

    // overriden rGetDaughterLocation() ? (not in AbstractCentreBasedDivisionRule.hpp)

    void OutputCellCentreBasedDivisionRuleInfo(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SproutingRule)

#endif /*SPROUTING_HPP_*/