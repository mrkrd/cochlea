all : _tw.so

bm_wave.o : bm_wave.c bm_wave.h
	gcc -fPIC -c bm_wave.c -o bm_wave.o

LCR4.o : LCR4.c LCR4.h
	gcc -fPIC -c LCR4.c -o LCR4.o

ihcrp.o : ihcrp.c ihcrp.h
	gcc -fPIC -c ihcrp.c -o ihcrp.o

# _bm.o : _bm.c
# 	gcc -fPIC `python-config --cflags` \
# 	-I `python -c "import numpy; print numpy.get_include()"` \
# 	-c _bm.c -o _bm.o

# _bm.so : _bm.o bm_wave.o LCR4.o ihcrp.o
# 	gcc -shared `python-config --ldflags` -o _bm.so \
# 	-L/usr/local/lib/pth \
# 	_bm.o bm_wave.o LCR4.o ihcrp.o

_tw.so : _tw.pyx bm_wave.o LCR4.o ihcrp.o
	cython _tw.pyx

	gcc -fPIC `python-config --cflags` \
	-I `python -c "import numpy; print numpy.get_include()"` \
	-c _tw.c

	gcc -shared `python-config --ldflags` \
	-L/usr/local/lib/pth \
	_tw.o bm_wave.o LCR4.o ihcrp.o \
	-o _tw.so

clean :
	rm -f _bm.so _bm.o bm_wave.o LCR4.o ihcrp.o
	rm -f _tw.so _tw.o _tw.c
	rm -f *.pyc
