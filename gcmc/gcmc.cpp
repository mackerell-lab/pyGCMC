/*
    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved
        Mingtian Zhao, Abhishek A. Kognole,
        Aoxiang Tao, Alexander D. MacKerell Jr.
    E-mail:
        zhaomt@outerbanks.umaryland.edu
        alex@outerbanks.umaryland.edu
*/


#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

extern "C" void vector_add_cuda(const float *a, const float *b, float *c, int n);
extern "C" void vector_sub_cuda(const float *a, const float *b, float *c, int n);


static PyObject *vector_add_py(PyObject *self, PyObject *args) {
    PyObject *py_a, *py_b, *py_c;
    if (!PyArg_ParseTuple(args, "OOO", &py_a, &py_b, &py_c)) {
        return NULL;
    }

    if (!PyArray_Check(py_a) || !PyArray_Check(py_b) || !PyArray_Check(py_c)) {
        PyErr_SetString(PyExc_TypeError, "Expected numpy arrays as input.");
        return NULL;
    }

    if (PyArray_TYPE((PyArrayObject *)py_a) != NPY_FLOAT32 ||
        PyArray_TYPE((PyArrayObject *)py_b) != NPY_FLOAT32 ||
        PyArray_TYPE((PyArrayObject *)py_c) != NPY_FLOAT32) {
        PyErr_SetString(PyExc_TypeError, "Expected float32 numpy arrays as input.");
        return NULL;
    }

    int ndim_a = PyArray_NDIM((PyArrayObject *)py_a);
    int ndim_b = PyArray_NDIM((PyArrayObject *)py_b);
    int ndim_c = PyArray_NDIM((PyArrayObject *)py_c);

    if (ndim_a != 1 || ndim_b != 1 || ndim_c != 1) {
        PyErr_SetString(PyExc_ValueError, "Expected 1-dimensional numpy arrays as input.");
        return NULL;
    }

    int n = (int)PyArray_DIM((PyArrayObject *)py_a, 0);

    if (PyArray_DIM((PyArrayObject *)py_b, 0) != n || PyArray_DIM((PyArrayObject *)py_c, 0) != n) {
        PyErr_SetString(PyExc_ValueError, "Input arrays should have the same length.");
        return NULL;
    }


    float *a = (float *)PyArray_DATA((PyArrayObject *)py_a);
    float *b = (float *)PyArray_DATA((PyArrayObject *)py_b);
    float *c = (float *)PyArray_DATA((PyArrayObject *)py_c);

    vector_add_cuda(a, b, c, n);

    Py_RETURN_NONE;
}


static PyObject *vector_sub_py(PyObject *self, PyObject *args) { // Add this new function
    PyObject *py_a, *py_b, *py_c;
    if (!PyArg_ParseTuple(args, "OOO", &py_a, &py_b, &py_c)) {
        return NULL;
    }

    if (!PyArray_Check(py_a) || !PyArray_Check(py_b) || !PyArray_Check(py_c)) {
        PyErr_SetString(PyExc_TypeError, "Expected numpy arrays as input.");
        return NULL;
    }

    if (PyArray_TYPE((PyArrayObject *)py_a) != NPY_FLOAT32 ||
        PyArray_TYPE((PyArrayObject *)py_b) != NPY_FLOAT32 ||
        PyArray_TYPE((PyArrayObject *)py_c) != NPY_FLOAT32) {
        PyErr_SetString(PyExc_TypeError, "Expected float32 numpy arrays as input.");
        return NULL;
    }

    int ndim_a = PyArray_NDIM((PyArrayObject *)py_a);
    int ndim_b = PyArray_NDIM((PyArrayObject *)py_b);
    int ndim_c = PyArray_NDIM((PyArrayObject *)py_c);

    if (ndim_a != 1 || ndim_b != 1 || ndim_c != 1) {
        PyErr_SetString(PyExc_ValueError, "Expected 1-dimensional numpy arrays as input.");
        return NULL;
    }

    int n = (int)PyArray_DIM((PyArrayObject *)py_a, 0);

    if (PyArray_DIM((PyArrayObject *)py_b, 0) != n || PyArray_DIM((PyArrayObject *)py_c, 0) != n) {
        PyErr_SetString(PyExc_ValueError, "Input arrays should have the same length.");
        return NULL;
    }

    float *a = (float *)PyArray_DATA((PyArrayObject *)py_a);
    float *b = (float *)PyArray_DATA((PyArrayObject *)py_b);
    float *c = (float *)PyArray_DATA((PyArrayObject *)py_c);

    vector_sub_cuda(a, b, c, n);

    Py_RETURN_NONE;
}


static PyMethodDef methods[] = {
    {"vector_add", vector_add_py, METH_VARARGS, "Add two vectors using CUDA"},
    {"vector_sub", vector_sub_py, METH_VARARGS, "Subtract two vectors using CUDA"}, 
    {NULL, NULL, 0, NULL}
};



PyMODINIT_FUNC PyInit_gpu(void) {
    PyObject *m;
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "gpu",          /* m_name */
        "Module doc",   /* m_doc */
        -1,             /* m_size */
        methods,        /* m_methods */
        NULL,           /* m_reload */
        NULL,           /* m_traverse */
        NULL,           /* m_clear */
        NULL,           /* m_free */
    };
    m = PyModule_Create(&moduledef);
    if (m == NULL)
        return NULL;
    import_array();  // Required for NumPy
    return m;
}
