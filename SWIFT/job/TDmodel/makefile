NESTLIBDIR = /Users/Danielyahalomi/Documents/Columbia/Research/Kipping/MultiNest/MultiNest_v3.12
LIBS = -L$(NESTLIBDIR) -llapack -lswift_gfortran -lnest3
FFLAGS = -O3 -DMPI
FC=mpif90

OBJFILES = params.o jasmine.o mandel.o plan.o ntrans_rmvs3.o transit.o like.o nestwrap.o main.o 

all: LRVplan 

%.o: %.f90
	$(FC) $(FFLAGS) -I$(NESTLIBDIR) -c $*.f90 
 
LRVplan: $(OBJFILES)
	$(FC) -I/usr/include/openmpi-x86_64/ -L./ -o ../tdmodel $(OBJFILES) \
	$(FFLAGS) $(LIBS) 

clean:
	rm -f *.o *.mod ../tdmodel
