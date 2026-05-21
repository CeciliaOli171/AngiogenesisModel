#include <../pybind11/include/pybind11/pybind11.h>
#include <../pybind11/include/pybind11/stl.h>
#include <vector>

namespace py = pybind11;

py::list ChasteDataForFetoflow(std::vector<double> v){
    py::list l;
    for(size_t i = 0; i < v.size(); ++i){
        l[i] = py::float_(v[i]);
    }
    return l;
}

std::vector<double> FetoflowDataForChaste(const py::iterable& l){
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