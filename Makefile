all : _pycat.so

_pycat.o : _pycat.pyx
	cython _pycat.pyx

	gcc -fPIC `python-config --cflags` \
	-I `python -c "import numpy; print numpy.get_include()"` \
	-c _pycat.c

catmodel_IHC.o : catmodel_IHC.c
	gcc -fPIC -c catmodel_IHC.c -o catmodel_IHC.o

catmodel_Synapse.o : catmodel_Synapse.c
	gcc -fPIC `python-config --cflags` \
	-c catmodel_Synapse.c -o catmodel_Synapse.o

complex.o : complex.c
	gcc -fPIC -c complex.c -o complex.o

_pycat.so : _pycat.o complex.o catmodel_Synapse.o catmodel_IHC.o
	gcc -shared `python-config --ldflags` \
	-o _pycat.so \
	_pycat.o catmodel_IHC.o complex.o catmodel_Synapse.o




clean :
	rm -f _pycat.so _pycat.o _pycat.h _pycat.c
	rm -f catmodel_IHC.o complex.o catmodel_Synapse.o
	rm -f *.pyc *.orig
