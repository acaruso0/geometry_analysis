#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "rdf.hpp"
#include "cdf.hpp"

PYBIND11_MODULE(corr_cpp, m) // name of the output module, variable
{
    m.doc() = "This is an example of C++ module called from python";
    m.def("couples", couples, "couples calculation"); // first one is the function in Python, second one is the function in C++
    m.def("rdf_dists", rdf_dists, "rdf distance calculation");
    m.def("triples", triples, "triples calculation");
    m.def("cdf_distang", cdf_distang, "cdf distance calculation");
}
