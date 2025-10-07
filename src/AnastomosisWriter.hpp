#ifndef ANASTOMOSISWRITER_HPP_
#define ANASTOMOSISWRITER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCellWriter.hpp"

/**
 * A class written using the visitor pattern for writing anastomosis events to file.
 *
 * The output file is called results.vizanastomosis by default. If VTK is switched on,
 * then the writer also specifies the VTK output for each cell, which is stored in
 * the VTK cell data "Anastomosis" by default.
 */

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AnastomosisWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{
private:
    /* serialization */
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

    /** 
     * Default constructor 
     */ 
    AnastomosisWriter();

    /** 
     * Overridden GetCellDataForVtkOutput() method.  
     *
     * Get a double equals to the type of anastomosis event associated with a cell.
     *
     * @param pCell a cell
     * @param pCellPopulation a pointer to the cell population owning the cell
     *
     * @return data associated with the cell
     */
    double GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /** 
     * Overridden VisitCell() method.  
     *
     * Visits a cell and write in a file the anastomosis event associated to it.
     *
     * @param pCell a cell
     * @param pCellPopulation a pointer to the cell population owning the cell
     */
    void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(AnastomosisWriter)

#endif /* ANASTOMOSISWRITER_HPP_ */
