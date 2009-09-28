all : _pycat.so

catmodel_IHC.o : catmodel_IHC.c
	gcc -fPIC -c catmodel_IHC.c -o catmodel_IHC.o

_pycat.so : _pycat.c catmodel_IHC.o
	gcc -fPIC `python-config --cflags` \
	-I `python -c "import numpy; print numpy.get_include()"` \
	-c _pycat.c -o _pycat.o

	gcc -shared `python-config --ldflags` `dsam-config --libs` -o _pycat.so \
	_pycat.o catmodel_IHC.o

clean :
	rm -f _pycat.so _pycat.o catmodel_IHC.o