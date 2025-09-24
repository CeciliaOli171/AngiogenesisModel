#ifndef BIRTHTIMECELLWRITER_HPP_
#define BIRTHTIMECELLWRITER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCellWriter.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class BirthTimeCellWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{
private:
    // Needed for serialization
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    // Constructor
    BirthTimeCellWriter();

    // Overridden GetCellDataForVtkOutput() method
    double GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    //Overridden VisitCell() method
    void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(BirthTimeCellWriter)

#endif /* BIRTHTIMECELLWRITER_HPP_ */
