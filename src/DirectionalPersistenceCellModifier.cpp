/*

Copyright (c) 2005-2016, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "DirectionalPersistenceCellModifier.hpp"
#include "Debug.hpp"

template<unsigned DIM>
DirectionalPersistenceCellModifier<DIM>::DirectionalPersistenceCellModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
    mDataLabels = std::vector<std::string>{ "fx", "fy", "fz" };
}

template<unsigned DIM>
DirectionalPersistenceCellModifier<DIM>::~DirectionalPersistenceCellModifier()
{
}

template<unsigned DIM>
void DirectionalPersistenceCellModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void DirectionalPersistenceCellModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    for ( typename AbstractCellPopulation<DIM>::Iterator cell_it = rCellPopulation.Begin(); 
                                                cell_it != rCellPopulation.End(); ++cell_it )
    {
        for (unsigned i = 0; i < DIM; i++)
        {
            (*cell_it)->GetCellData()->SetItem(mDataLabels[i], 0.0);
        }
        (*cell_it)->GetCellData()->SetItem("MovementAngle", 0.0);
    }
}

template<unsigned DIM>
void DirectionalPersistenceCellModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Check we are a node based cell population
    NodeBasedCellPopulation<DIM>* pPopulation = dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);
    // if ( (! (bool) pPopulation) | (DIM != 2) )
    // {
    //     EXCEPTION("Directional persistence cell modifier only implemented and checked for node based cell populations in 2D.\n");
    // }
    if ( ! (bool) pPopulation )
    {
        EXCEPTION("Directional persistence cell modifier only implemented and checked for node based cell populations in 2D.\n");
    }
    
    // Iterate over all the cells/nodes
    for ( typename AbstractCellPopulation<DIM>::Iterator cell_it = pPopulation->Begin(); 
                                                cell_it != pPopulation->End(); ++cell_it )
    {
        c_vector<double, DIM> force = pPopulation->GetNodeCorrespondingToCell(*cell_it)->rGetAppliedForce();
        double force_mag = norm_2(force);
        // force_mag = (force_mag == 0) ? 1.0:force_mag;
        for (unsigned i = 0; i < DIM; i++)
        {
            double unit_f = (force_mag < 1.0e-6) ? 0.0 : force[i]/force_mag;
            (*cell_it)->GetCellData()->SetItem(mDataLabels[i], unit_f);
        }
        double angle = atan2(force[DIM-1], force[0]);
        (*cell_it)->GetCellData()->SetItem("MovementAngle", angle);
    }
}

template<unsigned DIM>
void DirectionalPersistenceCellModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // Call parent method
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class DirectionalPersistenceCellModifier<1>;
template class DirectionalPersistenceCellModifier<2>;
template class DirectionalPersistenceCellModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DirectionalPersistenceCellModifier)

