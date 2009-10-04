#include "Python.h"
#include "numpy/arrayobject.h"

#include "catmodel_IHC.h"



static PyObject*
catmodel_IHC_wrap(PyObject* self, PyObject* args)
{

     PyObject *signal_arr, *signal_arg;
     double *signal_data;
     double cf, fs, cohc, cihc;

     PyObject *ihcout_arr;
     npy_intp dims[1];
     int signal_len;
     double *ihcout_data;

     int i;

     if (!PyArg_ParseTuple(args, "Odddd", \
     			   &signal_arg, &cf, &fs, \
     			   &cohc, &cihc))
     	  return NULL;


     /* Input: sound OR px */
     signal_arr = PyArray_FROM_OTF(signal_arg, NPY_DOUBLE, NPY_IN_ARRAY);
     signal_data = PyArray_DATA(signal_arr);
     signal_len = PyArray_DIM(signal_arr, 0);


     /* Output: ihcout OR ihcout_arr */
     dims[0] = signal_len;
     ihcout_arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
     ihcout_data = PyArray_DATA(ihcout_arr);
     /* TODO: use PyArray_Zeros instead */
     for (i=0; i<signal_len; i++) {
	  ihcout_data[i] = 0.0;
     }

     /* Run IHC function */
     SingleAN_IHC(signal_data, cf, 1, 1.0/fs, signal_len, cohc, cihc, ihcout_data);


     Py_DECREF(signal_arr);
     return ihcout_arr;
}



static PyMethodDef
PyCat_Methods[] =
{
     {"ihc", catmodel_IHC_wrap, METH_VARARGS, "IHC module."},
     {NULL, NULL, 0, NULL}
};



PyMODINIT_FUNC
init_pycat(void)
{
     (void)Py_InitModule("_pycat", PyCat_Methods);
     import_array();
}

