#include "Python.h"
#include "numpy/arrayobject.h"


double* ffGn(int N, double Hinput, double mu)
{

     double *data, *randNums;
     PyObject *modname, *mod, *mdict, *func, *result, *args;


     Py_Initialize();
     modname = PyString_FromString("ffGn");
     mod = PyImport_Import(modname);


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

     /* if (result) { */

     data = PyArray_DATA(result);
     randNums = malloc( N*sizeof(double) );

     memcpy(randNums, data, N*sizeof(double));


     Py_XDECREF(result);
     Py_XDECREF(args);
     /* Py_XDECREF(func); */
     /* Py_XDECREF(mdict); */
     Py_XDECREF(mod);
     Py_XDECREF(modname);

     Py_Finalize();

     return randNums;
}


double* pyResample(double *x, int orig_len, int new_len)
{
     PyObject *modname, *mod, *mdict, *func, *result, *args;
     PyObject *input_arr;
     double *input_data, *output_data;
     double *y;

     npy_intp dims[1];
     dims[0] = orig_len;


     Py_Initialize();
     modname = PyString_FromString("ffGn");
     mod = PyImport_Import(modname);
     printf("asdf\n");
     if (mod) {
	  printf("scipy loaded");
     }


     mdict = PyModule_GetDict(mod);
     func = PyDict_GetItemString(mdict, "resample"); /* borrowed reference */


     input_arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
     input_data = PyArray_DATA(input_arr);
     /* TODO: protect against buffer overflow */
     memcpy(input_data, x, orig_len*sizeof(double));

     args = PyTuple_New(3);
     PyTuple_SetItem(args, 0, input_arr);
     PyTuple_SetItem(args, 1, new_len);

     result = PyObject_CallObject(func, args);

     if (result) {
	  printf("*\n");
     }

     output_data = PyArray_DATA(result);
     y = malloc( new_len*sizeof(double) );

     memcpy(y, output_data, new_len*sizeof(double));


     Py_XDECREF(input_arr);
     Py_XDECREF(result);
     Py_XDECREF(args);
     /* Py_XDECREF(func); */
     Py_XDECREF(mdict);
     Py_XDECREF(mod);
     Py_XDECREF(modname);

     Py_Finalize();

     return y;
}
