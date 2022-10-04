#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <math.h>
#include <numpy/arrayobject.h>
#include <stdio.h>
#include <vector>

typedef npy_intp integer;
typedef npy_double real;

/*  wrapped cosine function */
static PyObject* davies_method(PyObject* self, PyObject* args)
{
    PyObject *arg1=NULL;
    PyArrayObject *lb = NULL, *nc = NULL;


    if (!PyArg_ParseTuple(args, "O", &arg1, &PyArray_Type))
        return NULL;

    lb = (PyArrayObject*) PyArray_FROM_OTF(arg1, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (lb == NULL)
        return NULL;

    integer nd = PyArray_NDIM(lb);

    if (nd != 1) {
        PyErr_SetString(PyExc_RuntimeError,
                        "The first argument needs to be a 1-dimensional vector");
        return NULL;
    }

    integer* dims = PyArray_SHAPE(lb);

    /*  construct the output from cos, from c double to python float */
    return Py_BuildValue("f", cos(nd));
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