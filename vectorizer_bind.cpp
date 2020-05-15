//
// Created by Sadeep on 27-Apr.
//

#include "vectorizer_bind.h"

using namespace std;

PyObject* freq_info_list_to_py(std::vector<freq_info>* list){
    PyObject* pylist = PyList_New(list->size());
    for(int c=0;c<list->size();c++){
        PyObject* t = PyTuple_New(3);
        PyTuple_SetItem(t, 0, PyLong_FromLong(list->at(c).start));
        PyTuple_SetItem(t, 1, PyLong_FromLong(list->at(c).end));
        PyObject* arr = PyList_New(32);
        for(int i=0;i<32;i++){
            PyList_SetItem(arr,i,PyLong_FromLong(list->at(c).freqs[i]));
        }
        PyTuple_SetItem(t, 2, arr);
        PyList_SetItem(pylist,c,t);
    }
    return pylist;

}
PyObject* freq_info_list_to_py4(std::vector<freq_info4>* list){
    PyObject* pylist = PyList_New(list->size());
    for(int c=0;c<list->size();c++){
        PyObject* t = PyTuple_New(3);
        PyTuple_SetItem(t, 0, PyLong_FromLong(list->at(c).start));
        PyTuple_SetItem(t, 1, PyLong_FromLong(list->at(c).end));
        PyObject* arr = PyList_New(128);
        for(int i=0;i<128;i++){
            PyList_SetItem(arr,i,PyLong_FromLong(list->at(c).freqs[i]));
        }
        PyTuple_SetItem(t, 2, arr);
        PyList_SetItem(pylist,c,t);
    }
    return pylist;

}

PyObject* _bind_split_n_count(PyObject *module,PyObject *args){

    char *raw_seq_internal;
    int n_internal;
    int length_mean_internal ;
    int length_stdev_internal;
    int count_internal;

    PyArg_ParseTuple(args, "siiii", &raw_seq_internal, &n_internal, &length_mean_internal,&length_stdev_internal,&count_internal);

    std::vector<freq_info> list;
    split_n_count(raw_seq_internal,n_internal,length_mean_internal,length_stdev_internal,count_internal,&list);

    return freq_info_list_to_py(&list);
}

PyObject* _bind_split_n_count_fasta(PyObject *module,PyObject *args){

    char *path;
    int length_mean_internal ;
    int length_stdev_internal;
    int count_internal;

    PyArg_ParseTuple(args, "siii", &path, &length_mean_internal, &length_stdev_internal, &count_internal);

    kseq_cpp::kseq_parser kseq(path);
    kseq.read_seq();

    std::vector<freq_info> list;
    split_n_count(kseq.seq.c_str(), kseq.seq.length(), length_mean_internal, length_stdev_internal, count_internal, &list);

    return freq_info_list_to_py(&list);
}
PyObject* _bind_split_n_count_fasta4(PyObject *module,PyObject *args){

    char *path;
    int length_mean_internal ;
    int length_stdev_internal;
    int count_internal;

    PyArg_ParseTuple(args, "siii", &path, &length_mean_internal, &length_stdev_internal, &count_internal);

    kseq_cpp::kseq_parser kseq(path);
    kseq.read_seq();

    std::vector<freq_info4> list;
    split_n_count4(kseq.seq.c_str(), kseq.seq.length(), length_mean_internal, length_stdev_internal, count_internal, &list);

    return freq_info_list_to_py4(&list);
}