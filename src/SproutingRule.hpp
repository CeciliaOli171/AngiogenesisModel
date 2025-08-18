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


// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCentreBasedCellPopulation;
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCentreBasedDivisionRule;

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class SproutingRule  : public AbstractCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>
{
friend class TestForcesModel;

private:

    double mMaxSproutingRate;
    double GetMaxSproutingRate();

    double mThresholdLength;

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
    SproutingRule(double maxSproutingRate = 0.08, double thresholdLength=2.0);

    // destructor 
    ~SproutingRule();

    // calculates the sprouting probability depending on the vegf concentration around the tip cell 
    virtual double GetSproutingProbability(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pParentCell);

    // calculates the set of all the vessel segment neighbours of the tip cell
    std::set<unsigned> GetVesselSegmentNeighbours(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, NodeBasedCellPopulation<SPACE_DIM>* p_node_population, CellPtr pParentCell); 

    // send the nth value of the neighbouring nodes set 
    unsigned GetNthNeighbourIndice(std::set<unsigned> neighbouring_node_indices, int n);
    
    // calculates the coordinates of the closest neighbour of a cell 
    std::pair<c_vector<double, SPACE_DIM>, unsigned> ClosestNeighbour(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pCell, std::set<unsigned> neighbouring_node_indices);

    bool IsBranchingCellNextToCell(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation,NodeBasedCellPopulation<SPACE_DIM>* p_node_population, CellPtr pParentCell);

    // calculates the coordinates of the perpendicular vector to the daughter position 
    c_vector<double, SPACE_DIM> PerpendicularDaughterDirection(c_vector<double, SPACE_DIM> daughter_direction, c_vector<double, SPACE_DIM> closest_neighbour);

    // calculates the daughter and parent position for a vessel segment division 
    // to modify and apply it to our model : for now no division of vessel segments 
    std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > CalculateCellDivisionVesselCellVector(CellPtr pParentCell, AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);

    // calculates the daughter and parent position for a tip cell division in the general case 
    std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > CalculateGrowthVector( AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, NodeBasedCellPopulation<SPACE_DIM>* p_node_population, CellPtr pParentCell);

    // calculates the daughter and parent position for a tip cell division 
    std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > CalculateCellDivisionTipCellVector(CellPtr pParentCell, AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);

    // overrides CalculateCellDivisionVector
    std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > CalculateCellDivisionVector(CellPtr pParentCell, AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);

    // overrides OutputCellCentreBasedDivisionRuleInfo
    void OutputCellCentreBasedDivisionRuleInfo(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SproutingRule)

#endif /*SPROUTINGRULE_HPP_*/