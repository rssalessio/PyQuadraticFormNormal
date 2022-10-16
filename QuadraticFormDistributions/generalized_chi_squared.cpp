#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "qfc.cpp"
#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

void check_error(int error_flag) {
    switch(error_flag) {
        case 1:
            throw std::range_error("Required accuracy not achieved");
        case 2:
            throw std::overflow_error("Round-off error possibly significant");
        case 3:
            throw std::invalid_argument("Invalid parameters");
        case 4:
            throw std::out_of_range("Unable to locate integration parameter");
        case 5:
            throw std::overflow_error("Out of memory");
        default:
            break;
    }
    return;
}

void check_parameters(
    py::array_t<double> x,
    py::array_t<double> coeff,
    py::array_t<double> nc,
    py::array_t<unsigned long int> df,
    double sigma,
    unsigned long int limit,
    double accuracy
) {
    if (coeff.size() != nc.size() || coeff.size() != df.size())
        throw std::length_error("coeff, nc and df need to have the same shape");
    else if (coeff.ndim() != 1 || x.ndim() != 1) 
        throw std::length_error("x, coeff, nc and df need to be 1-D list");
    else if (x.shape()[0] == 0)
        throw std::length_error("x cannot be empty");
    else if (accuracy <= 0)
        throw std::domain_error("Accuracy needs to be strictly positive");
    return;
}

std::tuple<py::array_t<double>, py::array_t<double>, int> davies_method(
    py::array_t<double> x,
    py::array_t<double> coeff,
    py::array_t<double> nc,
    py::array_t<unsigned long int> df,
    double sigma,
    unsigned long int limit,
    double accuracy) {

    check_parameters(x, coeff, nc, df, sigma, limit, accuracy);

    /* No pointer is passed, so NumPy will allocate the buffer */
    py::array_t<double> result = py::array_t<double>(x.shape()[0]);
    py::array_t<double> trace = py::array_t<double>(7);
    int ifault = 0;
    int r = coeff.shape()[0] + 1;

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

    check_error(ifault);

    return std::tuple(result, trace, ifault);
}

PYBIND11_MODULE(QuadraticFormDistributions, m) {
    m.doc() = R"pbdoc(
        QuadraticFormDistributions
        -----------------------
        .. currentmodule:: QuadraticFormDistributions
        .. autosummary::
           :toctree: _generate
           davies_method
    )pbdoc";

    m.def(
        "davies_method",
        &davies_method,
        R"pbdoc(
            Distribution function of quadratic forms in normal variables using Davies’s method.
            
            Parameters
            :param x: List of points to evaluate
            :param coeff: Coefficients of the chi^2 distributions.
            :param nc: Non-centrality parameters. Needs to be the same size as coeff.
            :param df: Degrees of freedom. Needs to be the same size as coeff.
            :param sigma: Std of the gaussian 
            :param limit: Maximum number of iterations
            :param accuracy: Desired accuracy

            Returns
            :result results: Probability of the evaluted points
            :result trace: Diagnostics
            :result fault: If 0, the algorithm has terminated successfully.        
        )pbdoc",
        py::arg("x"),
        py::arg("coeff"),
        py::arg("nc"),
        py::arg("df"),
        py::arg("sigma") = 1,
        py::arg("limit") = 10000,
        py::arg("accuracy") = 0.00001
        );

    #ifdef VERSION_INFO
        m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
    #else
        m.attr("__version__") = "dev";
    #endif
}