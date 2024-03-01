#include "CellPhaseCycle.hpp"

#include "CheckpointArchiveTypes.hpp"

#include "SmartPointers.hpp"
#include "Exception.hpp"

#include "AbstractCellCycleModel.hpp"
#include "AbstractSimpleGenerationalCellCycleModel.hpp"

#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"


CellPhaseCycle::CellPhaseCycle()
//   : AbstractCellCycleModel()
{
}

CellPhaseCycle::~CellPhaseCycle()
{
}

// overrides CreateCellCycleModel()
AbstractCellCycleModel* CellPhaseCycle::CreateCellCycleModel()
{
    CellPhaseCycle* p_model = new CellPhaseCycle();

    p_model->SetBirthTime(mBirthTime);
    p_model->SetMinimumGapDuration(mMinimumGapDuration);
    p_model->SetStemCellG1Duration(mStemCellG1Duration);
    p_model->SetGeneration(mGeneration);
    //p_model->GetG1Duration();
    //p_model->SetSDuration();
    //p_model->SetG2Duration();
    //p_model->SetMDuration();

    return p_model;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
//EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellPhaseCycle)
CHASTE_CLASS_EXPORT(CellPhaseCycle)