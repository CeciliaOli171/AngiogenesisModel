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

    double mOldNumNodes;
    double GetOldNumberofNodes();

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM,DIM> >(*this);
    }

public:

    DaughterCellModifier(double mOldNumNodes = 3);

    ~DaughterCellModifier();

    void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    double GetAngleFromVectors(c_vector<double,DIM> u, c_vector<double,DIM> v);

    double ClosestAngleVesselSegment(AbstractCellPopulation<DIM, DIM>& rCellPopulation, CellPtr cell_ptr, std::set<unsigned> neighbouring_node_indices);

    double OptimalAngleVesselSegment(AbstractCellPopulation<DIM, DIM>& rCellPopulation, CellPtr cell_ptr, std::set<unsigned> neighbouring_node_indices);

    double LengthVesselSegment(AbstractCellPopulation<DIM, DIM>& rCellPopulation, CellPtr cell_ptr, std::set<unsigned> neighbouring_node_indices);

    void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DaughterCellModifier)

#endif /*DAUGHTERCELLMODIFIER_HPP_*/
