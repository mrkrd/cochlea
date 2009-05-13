all : _bm.so

bm_wave.o : bm_wave.c bm_wave.h
	gcc -fPIC -c bm_wave.c -o bm_wave.o

LCR4.o : LCR4.c LCR4.h
	gcc -fPIC -c LCR4.c -o LCR4.o

_bm.o : _bm.c
	gcc -fPIC `python-config --cflags` \
	-I /usr/local/lib/python2.5/site-packages/numpy/core/include \
	-c _bm.c -o _bm.o


_bm.so : _bm.o bm_wave.o LCR4.o
	gcc -shared `python-config --ldflags` -o _bm.so _bm.o bm_wave.o LCR4.o

clean :
	rm -f _bm.so _bm.o bm_wave.o LCR4.o
