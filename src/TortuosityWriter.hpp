#ifndef TORTUOSITYWRITER_HPP_
#define TORTUOSITYWRITER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCellWriter.hpp"

/**
 * A class written using the visitor pattern for writing the tortuosity of the network to file.
 *
 * The output file is called results.viztortuosity by default. If VTK is switched on,
 * then the writer also specifies the VTK output for each cell, which is stored in
 * the VTK cell data "Tortuosity" by default.
 */

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class TortuosityWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
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
    TortuosityWriter();

    /** 
     * Overridden GetCellDataForVtkOutput() method.  
     *
     * Get a double equals to the tortuosity of the network.
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
     * Visits a cell and write in a file the tortuosity associated to it.
     *
     * @param pCell a cell
     * @param pCellPopulation a pointer to the cell population owning the cell
     */
    void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(TortuosityWriter)

#endif /* TORTUOSITYWRITER_HPP_ */
