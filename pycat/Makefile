all : _pycat_ihc.so

catmodel_IHC.o : catmodel_IHC.c
	gcc -fPIC -c catmodel_IHC.c -o catmodel_IHC.o

_pycat_ihc.so : _pycat_ihc.c catmodel_IHC.o
	gcc -fPIC `python-config --cflags` \
	-I `python -c "import numpy; print numpy.get_include()"` \
	-c _pycat_ihc.c -o _pycat_ihc.o

	gcc -shared `python-config --ldflags` `dsam-config --libs` -o _pycat_ihc.so \
	_pycat_ihc.o catmodel_IHC.o

clean :
	rm -f _pycat_ihc.so _pycat_ihc.o catmodel_IHC.o