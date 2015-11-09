cPROCESSOR := $(shell uname -m)

#ifeq ($(PROCESSOR),ia64)
  F90=gfortran
  FFLAGS=-g -C -O2 -ffree-form -I/opt/local/include -fbounds-check -ffpe-trap=invalid,overflow,zero #-fcheck=all
  FFLAGS2=$(FFLAGS)
  LDFLAGS=-L/opt/local/lib -lnetcdf -lnetcdff -framework vecLib

  .PHONY= 2d_test clean
#else
#  ifeq ($(PROCESSOR),x86_64)
#    F90=/usr/local/pgi/linux86-64/5.2/bin/pgf90
#    FFLAGS=-Kieee -fastsse #-Kieee #
#    LDFLAGS=-lnetcdf
#  else
#    F90=lf95
#    FFLAGS=-g --chk aesux --trap -I/usr/local/include #--o2  -I/usr/local/include #
#    LDFLAGS=-lnetcdf
#  endif
#endif

SOURCES= nDGmod.f90 \
				 coeff_update.f90 \
				 nDGsweep.f90 \

MODULES= testParameters.f90 \
         splitFunctions.f90 \
				 splitEvalFlux.f90 \
				 splitPositivityLimit.f90 \
				 unspFunctions.f90 \
				 unspPositivityLimit.f90 \
				 unspEvalFlux.f90 \
				
OBJECTS=$(SOURCES:.f90=.o)
MODOBJ=$(MODULES:.f90=.o)

SUBDIR = n4/

all: $(SOURCES) test_noddg_2d

2d_test: test_noddg_2d
	./test_noddg_2d

test_noddg_2d: $(MODOBJ) $(OBJECTS) unsplit_2d_nodal.f90
	$(F90) $(FFLAGS) $(MODOBJ) $(OBJECTS) unsplit_2d_nodal.f90 \
	         -o $@ $(LDFLAGS)

clean:
	rm -f *.o test_noddg_2d

%.o : %.f90
	$(F90) -c $(FFLAGS) $<

#	cp screen.out _ndgzhshu/n5/pRefine/rescale
#	cp screen.out _ndgsplzs/n5/pRefine/rescale
#	cp screen.out _ndgzhshu/$(SUBDIR)
#cp screen.out _ndgunlim/$(SUBDIR)
#cp screen.out _matrunc/$(SUBDIR)
#./test_noddg_2d 2>&1 | tee screen.out
