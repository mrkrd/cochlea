#include "Python.h"
#include "numpy/arrayobject.h"


double* ffGn(int N, double Hinput, double mu)
{

     double *data, *randNums;
     PyObject *mod, *mdict, *func, *result, *args;
     PyObject *output_arr;

     Py_Initialize();
     import_array();

     mod = PyImport_ImportModule("ffGn");

     /* if (mod) printf("success\n"); else printf("failure\n"); */

     mdict = PyModule_GetDict(mod);
     func = PyDict_GetItemString(mdict, "ffGn"); /* borrowed reference */

     /* if (func) { */
     /*      if (PyCallable_Check(func)) */
     /* args = PyTuple_New(3); */
     /* PyTuple_SetItem(args, 0, N); */
     /* PyTuple_SetItem(args, 1, Hinput); */
     /* PyTuple_SetItem(args, 2, mu); */

     args = Py_BuildValue( "(i,d,d)", N, Hinput, mu);

     result = PyObject_CallObject(func, args);

     output_arr = PyArray_FROM_OTF(result, NPY_DOUBLE, NPY_IN_ARRAY);

     data = PyArray_DATA(output_arr);
     randNums = malloc( N*sizeof(double) );
     memcpy(randNums, data, N*sizeof(double));


     Py_DECREF(output_arr);
     Py_XDECREF(result);
     Py_XDECREF(args);
     /* Py_XDECREF(func); */
     /* Py_XDECREF(mdict); */
     Py_XDECREF(mod);

     /* Py_Finalize(); */


     return randNums;
}


double* pyResample(double *x, int len, int p, int q)
{
     PyObject *modname, *mod, *mdict, *func, *result, *args;
     PyObject *input_arr;
     double *input_data;
     double *y;

     npy_intp dims[1];

     PyObject *output_arr;
     double *output_data;
     int new_len;

     int i;

     Py_Initialize();
     import_array();


     mod = PyImport_ImportModule("scipy.signal");
     if (!mod) {
     	  printf("Module not loaded.\n");
     	  exit(-1);
     }

     mdict = PyModule_GetDict(mod);
     if (!mdict) {
     	  printf("Dict not loaded.\n");
     	  exit(-1);
     }

     func = PyDict_GetItemString(mdict, "resample"); /* borrowed reference */
     if (!func) {
     	  printf("Function not found.\n");
     	  exit(-1);
     }



     dims[0] = len;
     input_arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
     input_data = PyArray_DATA(input_arr);
     /* TODO: protect against buffer overflow */
     memcpy(input_data, x, len*sizeof(double));



     new_len = (int)ceil(len * p / q);

     args = PyTuple_New(2);
     PyTuple_SetItem(args, 0, input_arr);
     PyTuple_SetItem(args, 1, PyInt_FromLong(new_len));


     result = PyObject_CallObject(func, args);

     output_arr = PyArray_FROM_OTF(result, NPY_DOUBLE, NPY_IN_ARRAY);
     output_data = PyArray_DATA(output_arr);

     y = malloc( new_len*sizeof(double) );
     memcpy(y, output_data, new_len*sizeof(double));


     Py_XDECREF(input_arr);
     Py_XDECREF(output_arr);
     Py_XDECREF(result);
     Py_XDECREF(args);
     /* Py_XDECREF(func); */
     /* Py_XDECREF(mdict); */
     Py_XDECREF(mod);

     /* Py_Finalize(); */

     /* for (i=0; i<new_len; i++) { */
     /* 	  y[i] = i; */
     /* 	  printf("%d: %f\n", i, y[i]); */
     /* } */

     /* printf("XXX %x\n", y); */

     return y;
}


double* pyRand(int len)
{
     PyObject *modname, *mod, *mdict, *func, *result, *args;
     double *y;

     PyObject *output_arr;
     double *output_data;

     Py_Initialize();
     import_array();

     mod = PyImport_ImportModule("numpy.random");

     mdict = PyModule_GetDict(mod);
     func = PyDict_GetItemString(mdict, "rand"); /* borrowed reference */


     args = PyTuple_New(1);
     PyTuple_SetItem(args, 0, PyInt_FromLong(len));

     result = PyObject_CallObject(func, args);
     output_arr = PyArray_FROM_OTF(result, NPY_DOUBLE, NPY_IN_ARRAY);

     output_data = PyArray_DATA(output_arr);
     y = malloc( len*sizeof(double) );
     memcpy(y, output_data, len*sizeof(double));


     Py_XDECREF(output_arr);
     Py_XDECREF(result);
     Py_XDECREF(args);
     /* Py_XDECREF(func); */
     /* Py_XDECREF(mdict); */
     Py_XDECREF(mod);

     /* Py_Finalize(); */


     return y;
}


