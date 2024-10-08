#include "RandomForce.hpp"

#include <random> // to generate random vector 

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"
#include "RandomNumberGenerator.hpp"
#include "Debug.hpp"
#include "SimulationTime.hpp"

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
    //TRACE("Begin Random Force");

    // initialisation 
    c_vector<double, DIM> randomforce; 
    c_vector<double, DIM> rand_vector; 

    double dt = SimulationTime::Instance()->GetTime(); 

    // we applied the force to the cell population corresponding 
    // here, the random force is applied to every cell 
    for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
         node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        // we collect the node data (index)
        unsigned node_index = (node_iter)->GetIndex();

        // Initialise a force vector
        for (unsigned i=0; i<DIM; i++)
        {
            double rand_norm = 1.0 - 2.0*(RandomNumberGenerator::Instance()->ranf());
            rand_vector[i] =  rand_norm;
        }

        double norm_rand_vector = norm_2(rand_vector);        
        if(norm_rand_vector != 0){
            randomforce = dt*(mSigma/norm_rand_vector)*rand_vector; // calculation of the force
        } else {
            randomforce = zero_vector<double>(DIM);
        }
        node_iter->AddAppliedForceContribution(randomforce);
        //PRINT_VECTOR(randomforce);
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
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RandomForce)