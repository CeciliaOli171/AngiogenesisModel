// #ifndef COUPLINGFETOFLOWCHASTE_HPP_
// #define COUPLINGFETOFLOWCHASTE_HPP_

// #include <pybind11/pybind11.h>
// #include <pybind11/stl.h>
// #include <vector>

// namespace py = pybind11;

// /**
//  * TO MODIFY!!!! 
//  * A class for converting c_vector into python list and vice verso to use the python library Fetoflow for flow calculations.
//  */

// template class<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// class CouplingFetoflowChaste : 
// {
// private:
//     // values ?

// public:

//     CouplingFetoflowChaste();

//     ~CouplingFetoflowChaste();

//     py::list ChasteDataForFetoflow(std::vector<double> v);

//     std::vector<double> FetoflowDataForChaste(const py::iterable& l);
    
// };

// #endif /* COUPLINGFETOFLOWCHASTE_ */