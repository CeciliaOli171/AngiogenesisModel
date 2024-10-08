#ifndef CONSECUTIVEBRANCHESWRITER_HPP_
#define CONSECUTIVEBRANCHESWRITER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCellWriter.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ConsecutiveBranchesWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{
private:
    // Needed for serialization
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
        archive & boost::serialization::base_object<AbstractCellWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    // Constructor
    ConsecutiveBranchesWriter();

    // Overridden GetCellDataForVtkOutput() method
    double GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    //Overridden VisitCell() method
    void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ConsecutiveBranchesWriter)

#endif /* CONSECUTIVEBRANCHESWRITER_HPP_ */
