#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <math.h>
#include <numpy/arrayobject.h>
#include <stdio.h>
#include <vector>

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

/*  wrapped cosine function */
static PyObject* davies_method(PyObject* self, PyObject* args)
{
    PyObject *arg1 = NULL, *arg2 = NULL, *arg3 = NULL;
    PyArrayObject *lb = NULL, *nc = NULL, *df = NULL;

    if (!PyArg_ParseTuple(args, "OOO", &arg1, &arg2, &arg3, &PyArray_Type))
        return NULL;

    lb = convert_double_array(arg1);
    nc = convert_double_array(arg2);
    df = convert_int_array(arg3);
    if (lb == NULL || nc == NULL || df == NULL)
        return NULL;

    npy_intp* lb_shape = PyArray_SHAPE(lb);
    npy_intp* nc_shape = PyArray_SHAPE(nc);
    npy_intp* df_shape = PyArray_SHAPE(df);

    if ((nc_shape[0] != lb_shape[0]) || (nc_shape[0] != df_shape[0])) {
        PyErr_SetString(PyExc_RuntimeError, "lb and nc and df do not have the same shape.");
        return NULL;
    }

    /*  construct the output from cos, from c double to python float */
    return Py_BuildValue("f", cos(lb_shape[0]));
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