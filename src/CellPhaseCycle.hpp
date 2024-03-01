#ifndef CELLPHASECYCLE_HPP_
#define CELLPHASECYCLE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "CheckpointArchiveTypes.hpp"

#include "RandomNumberGenerator.hpp"

#include "SmartPointers.hpp"
#include "Exception.hpp"

#include "AbstractCellCycleModel.hpp"
#include "AbstractSimpleGenerationalCellCycleModel.hpp"

#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

class CellPhaseCycle  : public AbstractSimpleGenerationalCellCycleModel
{

private:

    // allow to archive the force model object in a cell-based simulation 
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleGenerationalCellCycleModel>(*this);
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        archive & *p_gen;
        archive & p_gen;
    }

    void SetG1Duration()
    {
        assert(mpCell != NULL);

        double uniform_random_number_G1 = RandomNumberGenerator::Instance()->ranf();

        if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
        {
            mG1Duration = -log(uniform_random_number_G1)*GetStemCellG1Duration();
        }
        else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
        {
            mG1Duration = DBL_MAX;
        }
        else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
        {
            mG1Duration = DBL_MAX;
        }
        else
        {
            NEVER_REACHED;
        }
    }

    void SetSDuration()
    {
        assert(mpCell != NULL);

        double uniform_random_number_S = RandomNumberGenerator::Instance()->ranf();

        if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
        {
            mSDuration = -log(uniform_random_number_S)*GetSDuration();
        }
        else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
        {
            mSDuration = DBL_MAX;
        }
        else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
        {
            mSDuration = DBL_MAX;
        }
        else
        {
            NEVER_REACHED;
        }
    }

    void SetG2Duration()
    {
        assert(mpCell != NULL);

        double uniform_random_number_G2 = RandomNumberGenerator::Instance()->ranf();

        if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
        {
            mG2Duration = -log(uniform_random_number_G2)*GetG2Duration();
        }
        else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
        {
            mG2Duration = DBL_MAX;
        }
        else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
        {
            mG2Duration = DBL_MAX;
        }
        else
        {
            NEVER_REACHED;
        }
    }

    void SetMDuration()
    {
        assert(mpCell != NULL);

        double uniform_random_number_M = RandomNumberGenerator::Instance()->ranf();

        if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
        {
            mMDuration = -log(uniform_random_number_M)*GetMDuration();
        }
        else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
        {
            mMDuration = DBL_MAX;
        }
        else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
        {
            mMDuration = DBL_MAX;
        }
        else
        {
            NEVER_REACHED;
        }
    }

public:

    // constructor 
    CellPhaseCycle();

    // destructor 
    ~CellPhaseCycle();

    // overrides CreateCellCycleModel
    AbstractCellCycleModel* CreateCellCycleModel();

};

#include "SerializationExportWrapper.hpp"
//EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellPhaseCycle) 
CHASTE_CLASS_EXPORT(CellPhaseCycle)

#endif /*CELLPHASECYCLE_HPP_*/