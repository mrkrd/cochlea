all : _pycat.so

catmodel_IHC.o : catmodel_IHC.c
	gcc -fPIC -c catmodel_IHC.c -o catmodel_IHC.o

complex.o : complex.c
	gcc -fPIC -c complex.c -o complex.o

ffGn.o : ffGn.c
	gcc -fPIC `python-config --cflags` \
	-I `python -c "import numpy; print numpy.get_include()"` \
	-c ffGn.c -o ffGn.o


_pycat.so : _pycat.c catmodel_IHC.o complex.o ffGn.o
	gcc -fPIC `python-config --cflags` \
	-I `python -c "import numpy; print numpy.get_include()"` \
	-c _pycat.c -o _pycat.o

	gcc -shared `python-config --ldflags` `dsam-config --libs` -o _pycat.so \
	_pycat.o catmodel_IHC.o complex.o ffGn.o


ffGn_test : ffGn_test.c ffGn.o
	gcc `python-config --cflags` \
	-I `python -c "import numpy; print numpy.get_include()"` \
	-c ffGn_test.c -o ffGn_test.o

	gcc `python-config --ldflags` -o ffGn_test \
	ffGn_test.o ffGn.o

clean :
	rm -f _pycat.so _pycat.o catmodel_IHC.o complex.o ffGn.o
