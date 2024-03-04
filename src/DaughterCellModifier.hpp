#ifndef DAUGHTERCELLMODIFIER_HPP_
#define DAUGHTERCELLMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"

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
#include "DifferentiatedCellProliferativeType.hpp"
#include "AbstractCentreBasedDivisionRule.hpp"

template<unsigned DIM>
class DaughterCellModifier : public AbstractCellBasedSimulationModifier<DIM,DIM>
{
    double mOldNumNodes;

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM,DIM> >(*this);
    }

public:

    DaughterCellModifier(double mOldNumNodes = 3);

    ~DaughterCellModifier();

    double GetOldNumberofNodes();

    void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    double OptimalAngleVesselSegment(std::set<unsigned> neighbouring_node_indices, CellPtr cell_ptr, AbstractCellPopulation<DIM, DIM>& rCellPopulation);

    double LengthVesselSegment(std::set<unsigned> neighbouring_node_indices, CellPtr cell_ptr, AbstractCellPopulation<DIM, DIM>& rCellPopulation);

    void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DaughterCellModifier)

#endif /*DAUGHTERCELLMODIFIER_HPP_*/
