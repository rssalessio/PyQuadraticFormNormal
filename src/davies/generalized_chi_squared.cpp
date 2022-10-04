#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <math.h>
#include <numpy/arrayobject.h>
#include <stdio.h>
#include <vector>
#include "qfc.cpp"

typedef npy_intp integer;
typedef npy_double real;

static PyArrayObject* convert_double_array(PyObject *arg) {
    PyArrayObject* lb = (PyArrayObject*) PyArray_FROM_OTF(arg, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (lb == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Error while converting an array to an array of real numbers.");
        return NULL;
    }

    integer nd = PyArray_NDIM(lb);
    if (nd != 1) {
        PyErr_SetString(PyExc_RuntimeError, "The first argument needs to be a 1-dimensional vector");
        return NULL;
    }
    return lb;
}


static PyArrayObject* convert_int_array(PyObject *arg) {
    PyArrayObject* lb = (PyArrayObject*) PyArray_FROM_OTF(arg, NPY_INTP, NPY_ARRAY_IN_ARRAY);
    if (lb == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Error while converting an array to an array of real numbers.");
        return NULL;
    }

    integer nd = PyArray_NDIM(lb);
    if (nd != 1) {
        PyErr_SetString(PyExc_RuntimeError, "The first argument needs to be a 1-dimensional vector");
        return NULL;
    }
    return lb;
}

/*  check also https://github.com/cran/CompQuadForm/blob/master/src/qfc.cpp */
static PyObject* davies_method(PyObject* self, PyObject* args)
{
    /* Input arguments */
    PyObject *arg1 = NULL, *arg2 = NULL, *arg3 = NULL, *arg4 = NULL, *arg5 = NULL, *arg6 = NULL, *arg7 = NULL;
    PyArrayObject *coeff = NULL, *nc = NULL, *df = NULL, *x = NULL;

    /* Output arguments */
    PyObject *results = NULL, *trace = NULL;
    npy_intp ifault = 0;

    double sigma = 0, accuracy = 1e-4;
    npy_int limit = 10000;

    /*double *sigma, double *c1, int *lim1, double *acc,*/

    if (!PyArg_ParseTuple(args, "OOOOOOO", &arg1, &arg2, &arg3, &arg4, &arg5, &arg6, &arg7, &PyArray_Type))
        return NULL;

    coeff = convert_double_array(arg1);
    nc = convert_double_array(arg2);
    df = convert_int_array(arg3);
    x = convert_double_array(arg4);
    sigma = ((double*) arg5)[0];
    limit = ((npy_int*) arg6)[0];
    accuracy = ((double*) arg7)[0];
    

    if (coeff == NULL || nc == NULL || df == NULL || x == NULL || sigma == NULL || limit == NULL || accuracy == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Some of the parameters are null");
        return NULL;
    }

    npy_intp* coeff_shape = PyArray_SHAPE(coeff);
    npy_intp* nc_shape = PyArray_SHAPE(nc);
    npy_intp* df_shape = PyArray_SHAPE(df);
    npy_intp* x_shape = PyArray_SHAPE(x);

    if ((nc_shape[0] != coeff_shape[0]) || (nc_shape[0] != df_shape[0])) {
        PyErr_SetString(PyExc_RuntimeError, "lb and nc and df do not have the same shape.");
        return NULL;
    }

    if (x_shape[0] <= 0) {
        PyErr_SetString(PyExc_RuntimeError, "x cannot be empty");
        return NULL;
    }

    npy_intp r = coeff_shape[0];

    results = PyList_New(x_shape[0]);
    trace = PyList_New(7);

    qfc((double*) coeff, (double*) nc, (int *) df, (int *) &r, (double *) &sigma,  (double *) x, (int *) &limit, (double *) &accuracy, (double *) trace, (int *) &ifault, (double *) &results);

    /*  construct the output from cos, from c double to python float */
    return Py_BuildValue("OOif", results, trace, ifault, cos(1+accuracy));
}

/*  define functions in module */
static PyMethodDef methods[] =
{
     {"davies_method", davies_method, METH_VARARGS, "Distribution function of quadratic forms in normal variables using Davies’s method."},
     {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
    /* module initialization */
    /* Python version 3*/
    static struct PyModuleDef mod =
    {
        PyModuleDef_HEAD_INIT,
        "generalized_chi_squared", "Documentation",
        -1,
        methods
    };

    PyMODINIT_FUNC
    PyInit_generalized_chi_squared(void)
    {
        import_array();
        return PyModule_Create(&mod);
    }

#else

    /* module initialization */
    /* Python version 2 */
    PyMODINIT_FUNC
    initcos_module(void)
    {
        (void) Py_InitModule("generalized_chi_squared", methods);
    }

#endif