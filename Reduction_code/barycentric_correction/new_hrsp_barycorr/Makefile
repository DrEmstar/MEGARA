EXE = barycorr_HRSP

OBJ = angle.o astro.o config.o fits.o general.o numeric.o process.o vector.o
HEAD = angle.h astro.h config.h fits.h general.h numeric.h process.h vector.h

all: $(EXE)

angle.o: angle.c $(HEAD)
	gcc -Wall   -c angle.c
astro.o: astro.c $(HEAD)
	gcc -Wall   -c astro.c
config.o: config.c $(HEAD)
	gcc -Wall   -c config.c
fits.o: fits.c $(HEAD)
	gcc -Wall   -c fits.c
general.o: general.c $(HEAD)
	gcc -Wall   -c general.c
numeric.o: numeric.c $(HEAD)
	gcc -Wall   -c numeric.c
process.o: process.c $(HEAD)
	gcc -Wall   -c process.c
vector.o: vector.c $(HEAD)
	gcc -Wall   -c vector.c

barycorr_HRSP: barycorr_HRSP.o $(OBJ)
	gcc -Wall   -o barycorr_HRSP $(OBJ) barycorr_HRSP.o -lm
barycorr_HRSP.o: barycorr_HRSP.c $(HEAD)
	gcc -Wall   -c barycorr_HRSP.c

clean:
	/bin/rm -f $(EXE) *.o *~

