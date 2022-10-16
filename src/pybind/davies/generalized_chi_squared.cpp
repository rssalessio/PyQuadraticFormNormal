#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "qfc.cpp"
namespace py = pybind11;

std::tuple<py::array_t<double>, py::array_t<double>, int> davies_method(
        py::array_t<double> x,
        py::array_t<double> coeff,
        py::array_t<double> nc,
        py::array_t<unsigned long int> df,
        double sigma,
        unsigned long int limit,
        double accuracy) {

    if (coeff.size() != nc.size() || coeff.size() != df.size())
        throw std::length_error("coeff, nc and df need to have the same shape");
    else if (coeff.ndim() != 1 || x.ndim() != 1) 
        throw std::length_error("x, coeff, nc and df need to be 1-D list");
    else if (x.shape()[0] == 0)
        throw std::length_error("x cannot be empty");
    else if (accuracy <= 0)
        throw std::domain_error("Accuracy needs to be strictly positive");

    /* No pointer is passed, so NumPy will allocate the buffer */
    py::array_t<double> result = py::array_t<double>(x.shape()[0]);
    py::array_t<double> trace = py::array_t<double>(7);
    int ifault = 0;

    qfc(
        static_cast<double *>(coeff.request().ptr),
        static_cast<double *>(nc.request().ptr),
        static_cast<int *>(df.request().ptr),
        (int *) &r,
        &sigma,
        static_cast<double *>(x.request().ptr),
        (int* ) &limit,
        &accuracy,
        static_cast<double *>(trace.request().ptr),
        &ifault, 
        static_cast<double *>(result.request().ptr));


    return std::tuple(result, trace, ifault);
}

PYBIND11_MODULE(generalized_chi_squared, m) {
    m.doc() = "pybind11 generalized_chi_squared plugin"; // optional module docstring

    m.def(
        "davies_method",
        &davies_method,
        "Distribution function of quadratic forms in normal variables using Davies’s method.",
        py::arg("x"),
        py::arg("coeff"),
        py::arg("nc"),
        py::arg("df"),
        py::arg("sigma") = 1,
        py::arg("limit") = 10000,
        py::arg("accuracy") = 0.0001
        );
}