#include "Python.h"
#include "numpy/arrayobject.h"


double* ffGn(int N, double Hinput, double mu)
{

     double *data, *randNums;
     PyObject *modname, *mod, *mdict, *func, *result, *args;

     int i;


     Py_Initialize();
     modname = PyString_FromString("ffGn");
     mod = PyImport_Import(modname);



     if (mod) {
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
	  Py_XDECREF(mod);
	  Py_XDECREF(modname);

     }

     Py_Finalize();

     return randNums;
}
