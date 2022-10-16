#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <math.h>
#include <numpy/arrayobject.h>
#include <stdio.h>
#include <vector>
#include "qfc.cpp"


static PyArrayObject* convert_double_array(PyObject *arg) {
    PyArrayObject* lb = (PyArrayObject*) PyArray_FROM_OTF(arg, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (lb == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Error while converting an array to an array of real numbers.");
        return NULL;
    }

    npy_int nd = PyArray_NDIM(lb);
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

    npy_int nd = PyArray_NDIM(lb);
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
    PyObject *arg1 = NULL, *arg2 = NULL, *arg3 = NULL, *arg4 = NULL, *arg5 = NULL, *arg6 = NULL, *arg7 = NULL, *arg8 = NULL, *arg9 = NULL, *arg10 = NULL;
    PyArrayObject *coeff = NULL, *nc = NULL, *df = NULL, *x = NULL, *results = NULL, *trace = NULL, *ifault = NULL;


    // double sigma = 0, acc = 1e-4;
    // npy_int limit = 10000;

    /*double *sigma, double *c1, int *lim1, double *acc,*/

    if (!PyArg_ParseTuple(args, "OOOOOOOOOO", &arg1, &arg2, &arg3, &arg4, &arg5, &arg6, &arg7, &arg8, &arg9, &arg10, &PyArray_Type))
        return NULL;

    coeff = convert_double_array(arg1);
    nc = convert_double_array(arg2);
    df = convert_int_array(arg3);
    x = convert_double_array(arg4);
    results = convert_double_array(arg8);
    trace = convert_double_array(arg9);
    ifault = convert_int_array(arg10);

    double sigma = PyFloat_AsDouble(arg5);
    npy_int limit = PyLong_AsLong(arg6);
    double accuracy = PyFloat_AsDouble(arg7);

    if (coeff == NULL || nc == NULL || df == NULL || x == NULL || results == NULL || trace == NULL || ifault == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Some of the parameters are null");
        return NULL;
    }

    npy_intp* coeff_shape = PyArray_SHAPE(coeff);
    npy_intp* nc_shape = PyArray_SHAPE(nc);
    npy_intp* df_shape = PyArray_SHAPE(df);
    npy_intp* x_shape = PyArray_SHAPE(x);
    npy_intp r = coeff_shape[0];

    if ((nc_shape[0] != coeff_shape[0]) || (nc_shape[0] != df_shape[0])) {
        PyErr_SetString(PyExc_RuntimeError, "lb and nc and df do not have the same shape.");
        return NULL;
    } else if (x_shape[0] <= 0) {
        PyErr_SetString(PyExc_RuntimeError, "x cannot be empty");
        return NULL;
    }

    
    double _trace[6] ={1, 0, 0, 0, 0, 0};
    int _ifault = 0;
    double* _res = (double*) malloc(x_shape[0] * (sizeof(double)));

    qfc((double*) coeff, (double*) nc, (int *) df, (int *) &r, &sigma,  (double *) x, &limit, &accuracy, _trace, &_ifault, _res);

    for (int i=0; i < x_shape[0]; i++) {
        PyArray_SETITEM(results, (char*)PyArray_GETPTR1(results,i), PyFloat_FromDouble(_res[i]));
	}

    for (int i=0; i < 6; i++) {
        PyArray_SETITEM(trace, (char*)PyArray_GETPTR1(trace, i), PyFloat_FromDouble(_trace[i]));
    }

    PyArray_SETITEM(ifault, (char*)PyArray_GETPTR1(ifault, 0), PyFloat_FromDouble(_ifault));


    free(_res);

    /*  construct the output from cos, from c double to python float */
    return Py_BuildValue("OOO", ifault, trace, results);
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