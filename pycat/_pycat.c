#include "Python.h"
#include "numpy/arrayobject.h"

#include "catmodel_IHC.h"



static PyObject*
catmodel_IHC_wrap(PyObject* self, PyObject* args)
{

     PyObject *signal_arr, *signal_arg;
     double cf, nrep, binwidth, reptime, cohc, cihc;

     double *px;
     double totalstim;
     PyObject *ihcout_arr;
     npy_intp dims[1];
     double *ihcout;


     if (!PyArg_ParseTuple(args, "Odddddd", \
     			   &signal_arg, &cf, &nrep, &binwidth, \
     			   &reptime, &cohc, &cihc))
     	  return NULL;

     /* Input sound */
     signal_arr = PyArray_FROM_OTF(signal_arg, NPY_DOUBLE, NPY_IN_ARRAY);
     px = PyArray_DATA(signal_arr);

     /* Input: totalstim */
     totalstim = (npy_intp)floor((reptime*1e3)/(binwidth*1e3));

     /* Output: ihcout / ihcout_arr */
     dims[0] = totalstim;
     ihcout_arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
     ihcout = PyArray_DATA(ihcout_arr);

     /* Run IHC function */
     SingleAN(px, cf, nrep, binwidth, totalstim, cohc, cihc, ihcout);

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

