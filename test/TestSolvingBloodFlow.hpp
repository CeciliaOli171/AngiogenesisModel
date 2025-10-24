#ifndef TESTSOLVINGBLOODFLOW_HPP_
#define TESTSOLVINGBLOODFLOW_HPP_

#include <cxxtest/TestSuite.h>

#include "CellPropertyRegistry.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <iostream>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include "ProcessSpecificArchive.hpp"
#include "OutputFileHandler.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "ArchiveOpener.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "Debug.hpp"
#include "CheckpointArchiveTypes.hpp"

#include "AbstractCellBasedSimulationModifier.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "AbstractCellMutationState.hpp"
#include "AbstractCellPopulationBoundaryCondition.hpp"
#include "AbstractBoxDomainPdeModifier.hpp"
#include "AnastomosisWriter.hpp"
#include "AveragedSourceParabolicPde.hpp"
#include "Cell.hpp"
#include "CellBasedEventHandler.hpp"
#include "CellData.hpp"
#include "CellId.hpp"
#include "CellLabel.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FileComparison.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "ParabolicBoxDomainPdeModifier.hpp"
#include "ReplicatableVector.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "SmartPointers.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include "UniformSourceParabolicPde.hpp"

#include "ChemoForce.hpp"
#include "ChemoForceWithConstantVegf.hpp"
#include "ChemoForceWithAnalyticalApproximationPde.hpp"
#include "ChemoForceWithPdes.hpp"
#include "RandomForce.hpp"
#include "PersistenceForce.hpp"
#include "AngularForce.hpp"
#include "LinearMechanicalForceModified.hpp"
#include "DirectionalPersistenceCellModifier.hpp"

#include "SproutingRule.hpp"
#include "SproutingRuleWithConstantVegf.hpp"
#include "SproutingRuleWithAnalyticalApproximationPde.hpp"
#include "SproutingRuleWithPdes.hpp"
#include "DaughterCellModifier.hpp"
#include "PinnedCellsBoundaryCondition.hpp"

#include "BranchingSegmentMutationState.hpp"
#include "VesselTipMutationState.hpp"
#include "VesselSegmentMutationState.hpp"

#include "VegfEquationPde.hpp"
#include "VegfBoundaryCondition.hpp"
#include "MolecularConcentrationsBoxDomainPdeModifier.hpp"

#include "BloodFlow.hpp"
#include "BloodFlowAngiogenesis.hpp"

class TestSolvingBloodFlow : public CxxTest::TestSuite
{
};

#endif /*TESTSOLVINGBLOODFLOW_HPP_*/