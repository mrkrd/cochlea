#include "Python.h"
#include "numpy/arrayobject.h"

#include "bm_wave.h"
#include "LCR4.h"
#include "ihcrp.h"

static PyObject*
bm_init_wrap(PyObject* self, PyObject* args)
{
     double f_s;
     PyObject *Ls_arg, *Ls_arr;
     PyObject *Rs_arg, *Rs_arr;
     PyObject *Ct_arg, *Ct_arr;
     PyObject *Rbm_arg, *Rbm_arr;
     PyObject *Cbm_arg, *Cbm_arr;
     PyObject *Lbm_arg, *Lbm_arr;
     double Rh;
     double Lh;

     double *Ls_data;
     double *Rs_data;
     double *Ct_data;
     double *Rbm_data;
     double *Cbm_data;
     double *Lbm_data;


     if (!PyArg_ParseTuple(args, "dOOOOOOdd",
			   &f_s,
			   &Ls_arg,
			   &Rs_arg,
			   &Ct_arg,
			   &Rbm_arg,
			   &Cbm_arg,
			   &Lbm_arg,
			   &Rh,
			   &Lh))
     	  return NULL;


     /* Extract Numpy arrays from arguments. */
     Ls_arr = PyArray_FROM_OTF(Ls_arg, NPY_DOUBLE, NPY_IN_ARRAY);
     if (Ls_arr == NULL) return NULL;
     Rs_arr = PyArray_FROM_OTF(Rs_arg, NPY_DOUBLE, NPY_IN_ARRAY);
     if (Rs_arr == NULL) return NULL;
     Ct_arr = PyArray_FROM_OTF(Ct_arg, NPY_DOUBLE, NPY_IN_ARRAY);
     if (Ct_arr == NULL) return NULL;
     Rbm_arr = PyArray_FROM_OTF(Rbm_arg, NPY_DOUBLE, NPY_IN_ARRAY);
     if (Rbm_arr == NULL) return NULL;
     Cbm_arr = PyArray_FROM_OTF(Cbm_arg, NPY_DOUBLE, NPY_IN_ARRAY);
     if (Cbm_arr == NULL) return NULL;
     Lbm_arr = PyArray_FROM_OTF(Lbm_arg, NPY_DOUBLE, NPY_IN_ARRAY);
     if (Lbm_arr == NULL) return NULL;

     /* Get data pointers to Numpy arrays */
     Ls_data = PyArray_DATA(Ls_arr);
     Rs_data = PyArray_DATA(Rs_arr);
     Ct_data = PyArray_DATA(Ct_arr);
     Rbm_data = PyArray_DATA(Rbm_arr);
     Cbm_data = PyArray_DATA(Cbm_arr);
     Lbm_data = PyArray_DATA(Lbm_arr);

     /* Call C function */
     bm_init(f_s,
	     Ls_data,
	     Rs_data,
	     Ct_data,
	     Rbm_data,
	     Cbm_data,
	     Lbm_data,
	     Rh,
	     Lh);

     /* printf("in C: %f\n", Ls_data[0]); */

     /* TODO: goto fail */
     /* Memory leak possible */
     /* http://docs.scipy.org/doc/numpy/user/c-info.how-to-extend.html#example */
     Py_DECREF(Ls_arr);
     Py_DECREF(Rs_arr);
     Py_DECREF(Ct_arr);
     Py_DECREF(Rbm_arr);
     Py_DECREF(Cbm_arr);
     Py_DECREF(Lbm_arr);

     Py_RETURN_NONE;
}




static PyObject*
bm_wave_wrap(PyObject* self, PyObject* args)
{
     /* Input arrays */
     PyObject *input_arg, *input_arr;
     PyObject *ampl_corr_arg, *ampl_corr_arr;
     PyObject *Abm_arg, *Abm_arr;
     PyObject *Cbm_arg, *Cbm_arr;
     double *input_data;
     double *ampl_corr_data;
     double *Abm_data;
     double *Cbm_data;

     /* Output array */
     PyObject *xBM_arr;
     int nd = 2;
     npy_intp xBM_dims[2];
     double *xBM_data;

     int sample_num, section_num;
     int i;

     if (!PyArg_ParseTuple(args, "OOOO",
			   &input_arg,
			   &ampl_corr_arg,
			   &Abm_arg,
			   &Cbm_arg))
     	  return NULL;

     /* Convert input arguments into Numpy arrays */
     input_arr = PyArray_FROM_OTF(input_arg, NPY_DOUBLE, NPY_IN_ARRAY);
     if (input_arr == NULL) return NULL;
     ampl_corr_arr = PyArray_FROM_OTF(ampl_corr_arg, NPY_DOUBLE, NPY_IN_ARRAY);
     if (ampl_corr_arr == NULL) return NULL;
     Abm_arr = PyArray_FROM_OTF(Abm_arg, NPY_DOUBLE, NPY_IN_ARRAY);
     if (Abm_arr == NULL) return NULL;
     Cbm_arr = PyArray_FROM_OTF(Cbm_arg, NPY_DOUBLE, NPY_IN_ARRAY);
     if (Cbm_arr == NULL) return NULL;



     /* Collect data about dimensions */
     /* TODO: input_arr should be 1D squeezed array */
     sample_num = PyArray_DIMS(input_arr)[0];
     section_num = SECTIONS;


     /* Generate output Numpy array (xBM) */
     xBM_dims[0] = sample_num;
     xBM_dims[1] = section_num;
     xBM_arr = PyArray_SimpleNew(nd, xBM_dims, NPY_DOUBLE);

     /* Get data pointers */
     input_data = PyArray_DATA(input_arr);
     xBM_data = PyArray_DATA(xBM_arr);
     ampl_corr_data = PyArray_DATA(ampl_corr_arr);
     Abm_data = PyArray_DATA(Abm_arr);
     Cbm_data = PyArray_DATA(Cbm_arr);


     /* Simulation loop:  over time samples */
     for (i = 0; i < sample_num; i++) {
	  /* printf("i: %d\n", i); */
	  bm_wave(input_data[i],
	  	  &xBM_data[section_num*i],
	  	  ampl_corr_data,
	  	  Abm_data,
	  	  Cbm_data);
     }

     Py_DECREF(input_arr);
     Py_DECREF(ampl_corr_arr);
     Py_DECREF(Abm_arr);
     Py_DECREF(Cbm_arr);

     return xBM_arr;
}


static PyObject*
LCR4_init_wrap(PyObject* self, PyObject* args)
{
     double f_s;
     PyObject *freq_map_arg, *freq_map_arr;
     PyObject *Qmin_arg, *Qmin_arr;
     PyObject *SAT1_arg, *SAT1_arr;
     PyObject *SAT4_arg, *SAT4_arr;

     double *freq_map_data;
     double *Qmin_data;
     double *SAT1_data;
     double *SAT4_data;

     if (!PyArg_ParseTuple(args, "dOOOO",
			   &f_s,
			   &freq_map_arg,
			   &Qmin_arg,
			   &SAT1_arg,
			   &SAT4_arg))
     	  return NULL;


     /* Extract Numpy arrays from arguments. */
     freq_map_arr = PyArray_FROM_OTF(freq_map_arg, NPY_DOUBLE, NPY_IN_ARRAY);
     if (freq_map_arr == NULL) return NULL;
     Qmin_arr = PyArray_FROM_OTF(Qmin_arg, NPY_DOUBLE, NPY_IN_ARRAY);
     if (Qmin_arr == NULL) return NULL;
     SAT1_arr = PyArray_FROM_OTF(SAT1_arg, NPY_DOUBLE, NPY_IN_ARRAY);
     if (SAT1_arr == NULL) return NULL;
     SAT4_arr = PyArray_FROM_OTF(SAT4_arg, NPY_DOUBLE, NPY_IN_ARRAY);
     if (SAT4_arr == NULL) return NULL;


     /* Get data pointers to Numpy arrays */
     freq_map_data = PyArray_DATA(freq_map_arr);
     Qmin_data = PyArray_DATA(Qmin_arr);
     SAT1_data = PyArray_DATA(SAT1_arr);
     SAT4_data = PyArray_DATA(SAT4_arr);

     /* Call C function */
     LCR4_init(f_s,
	       freq_map_data,
	       Qmin_data,
	       SAT1_data,
	       SAT4_data);

     /* printf("in C: %f\n", Ls_data[0]); */

     /* TODO: goto fail */
     /* Memory leak possible */
     /* http://docs.scipy.org/doc/numpy/user/c-info.how-to-extend.html#example */
     Py_DECREF(freq_map_arr);
     Py_DECREF(Qmin_arr);
     Py_DECREF(SAT1_arr);
     Py_DECREF(SAT4_arr);

     Py_RETURN_NONE;
}


static PyObject*
LCR4_wrap(PyObject* self, PyObject* args)
{
     PyObject *xBM_arg, *xBM_arr;
     PyObject *Qmin_arg, *Qmin_arr;
     PyObject *Qmax_arg, *Qmax_arr;

     double *xBM_data;
     double *Qmin_data;
     double *Qmax_data;

     int sample_num, section_num;
     int i;

     if (!PyArg_ParseTuple(args, "OOO",
			   &xBM_arg,
			   &Qmin_arg,
			   &Qmax_arg))
     	  return NULL;


     /* Extract Numpy arrays from arguments. */
     xBM_arr = PyArray_FROM_OTF(xBM_arg, NPY_DOUBLE, NPY_OUT_ARRAY);
     if (xBM_arr == NULL) return NULL;
     Qmin_arr = PyArray_FROM_OTF(Qmin_arg, NPY_DOUBLE, NPY_IN_ARRAY);
     if (Qmin_arr == NULL) return NULL;
     Qmax_arr = PyArray_FROM_OTF(Qmax_arg, NPY_DOUBLE, NPY_IN_ARRAY);
     if (Qmax_arr == NULL) return NULL;


     /* Dimensions */
     sample_num = PyArray_DIMS(xBM_arr)[0];
     section_num = PyArray_DIMS(xBM_arr)[1];


     /* Get data pointers to Numpy arrays */
     xBM_data = PyArray_DATA(xBM_arr);
     Qmin_data = PyArray_DATA(Qmin_arr);
     Qmax_data = PyArray_DATA(Qmax_arr);


     /* Simulation loop:  over time samples */
     for (i = 0; i < sample_num; i++) {
	  /* printf("i: %d\n", i); */
	  LCR4(&xBM_data[section_num*i],
	       Qmin_data,
	       Qmax_data,
	       0,
	       section_num);
     }


     /* printf("in C: %f\n", Ls_data[0]); */

     /* TODO: goto fail */
     /* Memory leak possible */
     /* http://docs.scipy.org/doc/numpy/user/c-info.how-to-extend.html#example */
     Py_DECREF(Qmin_arr);
     Py_DECREF(Qmax_arr);

     return xBM_arr;
}



static PyObject*
ihcrp_init_wrap(PyObject* self, PyObject* args)
{
     double f_s;

     if (!PyArg_ParseTuple(args, "d",
			   &f_s))
     	  return NULL;


     /* Call C function */
     ihcrp_init(f_s);

     Py_RETURN_NONE;
}

static PyObject*
ihcrp_wrap(PyObject* self, PyObject* args)
{
     /* Input arrays */
     PyObject *xBM_arg, *xBM_arr;
     PyObject *ciliaCouplingGain_arg, *ciliaCouplingGain_arr;
     double *xBM_data;
     double *ciliaCouplingGain_data;


     /* Output array */
     PyObject *uIHC_arr;
     int nd = 2;
     npy_intp uIHC_dims[2];
     double *uIHC_data;

     int sample_num, section_num;
     int i;

     if (!PyArg_ParseTuple(args, "OO",
			   &xBM_arg,
			   &ciliaCouplingGain_arg))
     	  return NULL;

     /* Convert input arguments into Numpy arrays */
     xBM_arr = PyArray_FROM_OTF(xBM_arg, NPY_DOUBLE, NPY_IN_ARRAY);
     if (xBM_arr == NULL) return NULL;
     ciliaCouplingGain_arr = PyArray_FROM_OTF(ciliaCouplingGain_arg, NPY_DOUBLE, NPY_IN_ARRAY);
     if (ciliaCouplingGain_arr == NULL) return NULL;


     /* Collect data about dimensions */
     sample_num = PyArray_DIMS(xBM_arr)[0];
     section_num = PyArray_DIMS(xBM_arr)[1];


     /* Generate output Numpy array (uIHC) */
     uIHC_dims[0] = sample_num;
     uIHC_dims[1] = section_num;
     uIHC_arr = PyArray_SimpleNew(nd, uIHC_dims, NPY_DOUBLE);

     /* Get data pointers */
     xBM_data = PyArray_DATA(xBM_arr);
     uIHC_data = PyArray_DATA(uIHC_arr);
     ciliaCouplingGain_data = PyArray_DATA(ciliaCouplingGain_arr);


     /* Simulation loop:  over time samples */
     for (i = 0; i < sample_num; i++) {
	  /* printf("i: %d\n", i); */
	  ihcrp(&uIHC_data[section_num*i],
		&xBM_data[section_num*i],
		ciliaCouplingGain_data);
     }

     Py_DECREF(xBM_arr);
     Py_DECREF(ciliaCouplingGain_arr);

     return uIHC_arr;
}




static PyMethodDef
BM_Methods[] =
{
     {"bm_init", bm_init_wrap, METH_VARARGS, "Call before bm_wave()."},
     {"bm_wave", bm_wave_wrap, METH_VARARGS, "BM processing."},
     {"LCR4_init", LCR4_init_wrap, METH_VARARGS, "Call before LCR4_wrap()."},
     {"LCR4", LCR4_wrap, METH_VARARGS, "LCR4 module."},
     {"ihcrp_init", ihcrp_init_wrap, METH_VARARGS, "Call before ihcrp_wrap()."},
     {"ihcrp", ihcrp_wrap, METH_VARARGS, "ihcrp module."},
     {NULL, NULL, 0, NULL}
};



PyMODINIT_FUNC
init_bm(void)
{
     (void)Py_InitModule("_bm", BM_Methods);
     import_array();
}

