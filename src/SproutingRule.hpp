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

private:

    double mPsprout;
    double mDirectionCoeff;
    double mPerpendicularDirectionCoeff;
    double GetSproutingProbability();
    double GetDirectionCoefficient();
    double GetPerpendicularDirectionCoefficient();

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
    SproutingRule(double Psprout = 0.1, double directioncoeff = 0.1, double perpendiculardirectioncoeff = 0.1);

    // destructor 
    ~SproutingRule();

    // send the nth value of the neighbouring nodes set 
    unsigned GetNthNeighbourIndice(std::set<unsigned> neighbouring_node_indices, int n);

    // calculates the coordinates of a random neighbour of a cell 
    std::tuple<CellPtr, unsigned, c_vector<double, SPACE_DIM>> RandomNeighbour(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pCell, std::set<unsigned> neighbouring_node_indices);
    
    // calculates the coordinates of the closest neighbour of a cell 
    c_vector<double, SPACE_DIM> ClosestNeighbour(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pCell, std::set<unsigned> neighbouring_node_indices);

    // calculates the vector between the two closest neighbours of a cell 
    c_vector<double, SPACE_DIM> TwoClosestNeighbours(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, CellPtr pCell, std::set<unsigned> neighbouring_node_indices);

    // overrides CalculateCellDivisionVector
    std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > CalculateCellDivisionVector(CellPtr pParentCell, AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);

    // overrides OutputCellCentreBasedDivisionRuleInfo
    void OutputCellCentreBasedDivisionRuleInfo(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SproutingRule)

#endif /*SPROUTINGRULE_HPP_*/