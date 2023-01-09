FORTRAN=gfortran
#COPTS= -c  -O -g -pedantic -fbounds-check 
COPTS= -c  -O3 -g 
#-c -checkbound -O -r8 -v
PROG=kalman
SUFFIXES=.f90.o
.SUFFIXES : .f90 .o
$(SUFFIXES):
	$(FORTRAN) $(COPTS) $*.f90

OBJS  = longr.o \
init.o\
fsource.o \
kalman.o 

export LDFLAGS="-L/opt/homebrew/opt/lapack/lib"
export CPPFLAGS="-I/opt/homebrew/opt/lapack/include"

$(PROG): $(OBJS)
	$(FORTRAN) $(OBJS) -llapack -lblas -o $(PROG)

clean:;	rm -f *.o *.mod objets/*



