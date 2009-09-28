all : _pycat_ihc.so

_pycat_ihc.so : _pycat_ihc.c
	gcc -fPIC `python-config --cflags` \
	-I `python -c "import numpy; print numpy.get_include()"` \
	-c _pycat_ihc.c -o _pycat_ihc.o

	gcc -shared `python-config --ldflags` `dsam-config --libs` -o _pycat_ihc.so _pycat_ihc.o

clean :
	rm -f _pycat_ihc.so _pycat_ihc.o