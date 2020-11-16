CC=armclang
F90=armflang

CFLAGS= -fopenmp -mcpu=native -O3 -ffp-contract=fast
FFLAGS= -fopenmp -mcpu=native -O3 -ffp-contract=fast

all: test-open.x test-arm.x

test-open.x: main.o routines.o
        $(F90) $(FFLAGS)  -o test-open.x main.o routines.o  ${HOME}/OpenBLAS/install/arm-20.0/0.3.9/lib/libopenblas.a -ldl -lm

test-arm.x: main.o routines.o
        $(F90) $(FFLAGS) -o test-arm.x main.o routines.o   ${BLAS_STATIC} ${LAPACK_STATIC}  -ldl -lm

clean:
	\rm *.o *.x *.mod

%.o: %.c
        $(CC) $(CFLAGS) -g -c $< -o $@

%.o: %.f90
        $(F90) $(FFLAGS) -g -c $< -o $@

