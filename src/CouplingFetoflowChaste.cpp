#include "CouplingFetoflowChaste.hpp"
#include "ReplicatableVector.hpp"

#include "Warnings.hpp"
#include "Debug.hpp"

#include <pybind11/pybind.h>
#include <pybind11/stl.h>
#include <vector>

namespace py = pybind11;

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CouplingFetoflowChaste<ELEMENT_DIM, SPACE_DIM>::CouplingFetoflowChaste()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CouplingFetoflowChaste<ELEMENT_DIM, SPACE_DIM>::~CouplingFetoflowChaste()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
py::list CouplingFetoflowChaste<ELEMENT_DIM, SPACE_DIM>::ChasteDataForFetoflow(std::vector<double> v){
    py::list l;
    l.resize(v.size());
    for(size_t i = 0; i < v.size(); ++i){
        l[i] = py::float_(v[i]);
    }
    return l;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> CouplingFetoflowChaste<ELEMENT_DIM, SPACE_DIM>::FetoflowDataForChaste(const py::iterable& l){
    std::vector<double> v;
    v.reserve(py::len(l));
    for (auto item : l){
        v.push_back(py::cast<double>(item));
    }
    return v;
}

// Modules 
PYBIND11_MODULE(CouplingFetoflowChaste, m){
    m.def("ChasteDataForFetoflow", &ChasteDataForFetoflow, py::arg("v"), "Convert a c++ vector into a python list");
    m.def("FetoflowDataForChaste", &FetoflowDataForChaste, py::arg("iterable"), "Convert a python list into a c++ vector");
}

// Explicit instantiation
template class CouplingFetoflowChaste<1,1>;
template class CouplingFetoflowChaste<1,2>;
template class CouplingFetoflowChaste<2,2>;
template class CouplingFetoflowChaste<1,3>;
template class CouplingFetoflowChaste<2,3>;
template class CouplingFetoflowChaste<3,3>;