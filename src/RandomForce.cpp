#include "RandomForce.hpp"

#include <random> // to generate random vector 

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"

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
    // we create the vector force 
    c_vector<double, DIM> randomforce = zero_vector<double>(DIM); 
    c_vector<double, DIM> rand_vector = zero_vector<double>(DIM); // initialization of a random vector 

    // we applied the force to the cell population corresponding 
    // here, the random force is applied to every cell 
    unsigned node_index = 0;
    for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
         node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        CellPtr pCell = rCellPopulation.GetCellUsingLocationIndex(node_index);
        if (pCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>()){
            //if(DIM == 2){
            //    rand_vector(0) =(float) std::rand()/RAND_MAX;
            //    rand_vector(1) =(float) std::rand()/RAND_MAX;
            // } else if(DIM == 3){
            //     rand_vector(0) =(float) std::rand()/RAND_MAX;
            //     rand_vector(1) =(float) std::rand()/RAND_MAX;
            //     rand_vector(2) =(float) std::rand()/RAND_MAX;
            // } else{
            //     rand_vector(0) =(float) std::rand()/RAND_MAX;
            // }

            // Initialise a force vector
            c_vector<double, DIM> rand_vector;
            for (unsigned i=0; i<DIM; i++)
            {
                /*
                    * The force on this cell is scaled with the timestep such that when it is
                    * used in the discretised equation of motion for the cell, we obtain the
                    * correct formula
                    *
                    * x_new = x_old + sqrt(2*D*dt)*W
                    *
                    * where W is a standard normal random variable.
                    * Note this is assuming a drag coefficient of 1.0.
                    */
                double rand_norm = 1.0 - 2.0*RandomNumberGenerator::Instance()->ranf();
                rand_vector[i] =  rand_norm;
            }

            double norm_rand_vector = norm_2(rand_vector);        
            if(norm_rand_vector != 0){
                randomforce = (mSigma/norm_rand_vector)*rand_vector; // calculation of the force
            }    
            node_iter->AddAppliedForceContribution(randomforce);
            }
        ++node_index;
    }
}

template<unsigned DIM>
void RandomForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    double sigma = GetRandomSensitivity();
    *rParamsFile << "\t\t\t<Sensitivity to Random Fluctuations>" << sigma << "</Sigma>\n";

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
//CHASTE_CLASS_EXPORT(RandomForce)