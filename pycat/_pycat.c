#include "Python.h"
#include "numpy/arrayobject.h"

#include "catmodel_IHC.h"

static PyMethodDef
DSAM_Methods[] =
{
     {"pycat_ihc", catmodel_IHC_wrap, METH_VARARGS, "IHC module."},
     {NULL, NULL, 0, NULL}
};



PyMODINIT_FUNC
init_dsam(void)
{
     (void)Py_InitModule("_dsam", DSAM_Methods);
     import_array();
}


