all : _catmodel.so

_pycat.o: _pycat.pyx
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

ffGn.o : ffGn.c
	gcc -fPIC `python-config --cflags` \
	-I `python -c "import numpy; print numpy.get_include()"` \
	-c ffGn.c -o ffGn.o


_catmodel.so : _catmodel.c catmodel_IHC.o complex.o ffGn.o  _pycat.o catmodel_Synapse.o
	gcc -fPIC `python-config --cflags` \
	-I `python -c "import numpy; print numpy.get_include()"` \
	-c _catmodel.c -o _catmodel.o

	gcc -shared `python-config --ldflags` -o _catmodel.so \
	_catmodel.o catmodel_IHC.o complex.o _pycat.o ffGn.o catmodel_Synapse.o \
	-L/usr/local/lib/pth



ffGn_test : ffGn_test.c ffGn.o
	gcc `python-config --cflags` \
	-I `python -c "import numpy; print numpy.get_include()"` \
	-c ffGn_test.c -o ffGn_test.o

	gcc `python-config --ldflags` -o ffGn_test \
	ffGn_test.o ffGn.o

clean :
	rm -f _catmodel.so _catmodel.o catmodel_IHC.o complex.o ffGn.o
	rm -f _pycat.o _pycat.c _pycat.h catmodel_Synapse.o
