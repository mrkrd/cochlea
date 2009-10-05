#include "Python.h"
#include "numpy/arrayobject.h"


double* ffGn(int N, double Hinput, double mu)
{

     double *data, *randNums;
     PyObject *modname, *mod, *mdict, *func, *result, *args;

     Py_Initialize();

     mod = PyImport_ImportModule("ffGn");

     mdict = PyModule_GetDict(mod);
     func = PyDict_GetItemString(mdict, "ffGn"); /* borrowed reference */

     /* if (func) { */
     /*      if (PyCallable_Check(func)) */
     /* args = PyTuple_New(3); */
     /* PyTuple_SetItem(args, 0, N); */
     /* PyTuple_SetItem(args, 1, Hinput); */
     /* PyTuple_SetItem(args, 2, mu); */

     args = Py_BuildValue( "(i,f,f)", N, Hinput, mu);

     result = PyObject_CallObject(func, args);

     data = PyArray_DATA(result);
     randNums = malloc( N*sizeof(double) );
     memcpy(randNums, data, N*sizeof(double));


     Py_XDECREF(result);
     Py_XDECREF(args);
     Py_XDECREF(func);
     Py_XDECREF(mdict);
     Py_XDECREF(mod);

     /* Py_Finalize(); */


     return randNums;
}


double* pyResample(double *x, int len, int p, int q)
{
     PyObject *modname, *mod, *mdict, *func, *result, *args;
     PyObject *input_arr;
     double *input_data, *result_data;
     double *y;

     npy_intp dims[1];

     PyObject *output_arr;
     double *output_data;
     int new_len;

     Py_Initialize();
     import_array();

     mod = PyImport_ImportModule("scipy.signal");

     mdict = PyModule_GetDict(mod);
     func = PyDict_GetItemString(mdict, "resample"); /* borrowed reference */



     dims[0] = len;
     input_arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
     input_data = PyArray_DATA(input_arr);
     /* TODO: protect against buffer overflow */
     memcpy(input_data, x, len*sizeof(double));



     new_len = ceil(len * p / q);
     printf("%d %d\n", len, new_len);
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
     Py_XDECREF(func);
     Py_XDECREF(mdict);
     Py_XDECREF(mod);
     Py_XDECREF(modname);

     /* Py_Finalize(); */


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
     Py_XDECREF(func);
     Py_XDECREF(mdict);
     Py_XDECREF(mod);
     Py_XDECREF(modname);

     /* Py_Finalize(); */


     return y;
}


