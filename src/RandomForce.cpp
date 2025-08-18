#include "RandomForce.hpp"

#include <random> // to generate random vector 

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "RandomNumberGenerator.hpp"
#include "Debug.hpp"
#include "SimulationTime.hpp"
#include "BranchingCellMutationState.hpp"

template<unsigned DIM>
RandomForce<DIM>::RandomForce(double sigma)
    : AbstractForce<DIM>()
{
    assert(sigma > 0);
    mSigma = sigma;
}

template<unsigned DIM>
RandomForce<DIM>::~RandomForce()
{
}

template<unsigned DIM>
double RandomForce<DIM>::GetRandomSensitivity()
{ 
    return mSigma;
}

// overrides AddForceContribution()
template<unsigned DIM>
void RandomForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // TRACE("Begin Random Force");

    // we applied the force to the cell population corresponding 
    // here, the random force is applied to every cell 
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // initialisation 
        c_vector<double, DIM> randomforce = zero_vector<double>(DIM); 
        c_vector<double, DIM> rand_vector; 

        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

        // Initialise a force vector
        for (unsigned i=0; i<DIM; i++)
        {
            double rand_norm = 1.0 - 2.0*(RandomNumberGenerator::Instance()->ranf());
            rand_vector[i] =  rand_norm;
        }

        double norm_rand_vector = norm_2(rand_vector);        
        if(norm_rand_vector != 0){
            randomforce = (mSigma/norm_rand_vector)*rand_vector; // calculation of the force
        }

        rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(randomforce);
        
        // test time force 
        //for(int i = 0; i < 100; ++i){rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(randomforce);}
    }

    //TRACE("End Random Force");
}

template<unsigned DIM>
void RandomForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<Sensitivity to Random Fluctuations>" << mSigma << "</Sigma>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class RandomForce<1>;
template class RandomForce<2>;
template class RandomForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RandomForce);