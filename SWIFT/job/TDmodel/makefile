NESTLIBDIR = /Users/Danielyahalomi/Documents/Columbia/Research/Kipping/MultiNest/MultiNest_v2.17/
LIBS = -L/Users/Danielyahalomi/Documents/Columbia/Research/Kipping/SWIFT/ -lswift -L$(NESTLIBDIR) -lnest3
FFLAGS = -O2
FC=gfortran

OBJFILES = params.o jasmine.o mandel.o plan.o ntrans_rmvs3.o transit.o like.o nestwrap.o main.o 

all: LRVplan 

%.o: %.f90
	$(FC) $(FFLAGS) -I$(NESTLIBDIR) -c $*.f90  
 
LRVplan: $(OBJFILES)
	$(FC) -o ../tdmodel.x $(OBJFILES) \
	$(FFLAGS) $(LIBS) -llapack  

clean:
	rm -f *.o *.mod ../tdmodel.x
