//
// Created by Sadeep on 27-Apr.
//

#ifndef PRE_PROCESSOR_VECTORIZER_BIND_H
#define PRE_PROCESSOR_VECTORIZER_BIND_H

#include <Python.h>
#include "kseq_cpp.h"
#include "vectorizer.h"


PyObject* tanh_impl(PyObject*, PyObject* o) {
    char* str1,* str2;
    if (!PyArg_ParseTuple(o,"ss", &str1, &str2))
    {
        return nullptr;
    }
    std::cout << str1 << std::endl << str2 << std::endl;
    return Py_None;
}

PyObject* _bind_split_n_count(PyObject *module,PyObject *args);
PyObject* _bind_split_n_count_fasta(PyObject *module,PyObject *args);;
PyObject* _bind_split_n_count_fasta4(PyObject *module,PyObject *args);

static PyMethodDef methods[] = {

        { "split_n_count", (PyCFunction)_bind_split_n_count, METH_VARARGS, nullptr },
        { "split_n_count_fasta", (PyCFunction)_bind_split_n_count_fasta, METH_VARARGS, nullptr },
        { "split_n_count_fasta4", (PyCFunction)_bind_split_n_count_fasta4, METH_VARARGS, nullptr },
        { nullptr, nullptr, 0, nullptr }
};

static PyModuleDef module = {
        PyModuleDef_HEAD_INIT,
        "vectorizer",
        "Python Package for Genome Vectorizer",
        0,
        methods
};

PyMODINIT_FUNC PyInit_vectorizer() {
    return PyModule_Create(&module);
}

#endif //PRE_PROCESSOR_VECTORIZER_BIND_H
