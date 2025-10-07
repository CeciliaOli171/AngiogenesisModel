#ifndef CONSECUTIVEBRANCHESWRITER_HPP_
#define CONSECUTIVEBRANCHESWRITER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCellWriter.hpp"

/**
 * A class written using the visitor pattern for writing number of consecutive branches to file.
 *
 * The output file is called results.vizconsecutivebranches by default. If VTK is switched on,
 * then the writer also specifies the VTK output for each cell, which is stored in
 * the VTK cell data "ConsecutiveBranches" by default.
 */

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ConsecutiveBranchesWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
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
    ConsecutiveBranchesWriter();

    /** 
     * Overridden GetCellDataForVtkOutput() method.  
     *
     * Get a double equals to the number of consecutive branches leading to a cell.
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
     * Visits a cell and write in a file the number of consecutive branches associated to it.
     *
     * @param pCell a cell
     * @param pCellPopulation a pointer to the cell population owning the cell
     */
    void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ConsecutiveBranchesWriter)

#endif /* CONSECUTIVEBRANCHESWRITER_HPP_ */
