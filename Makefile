# Yippee!  This is my very own makefile.
# Jamison A. Smith, AKA JAS
# April 26, 2007

FORTRAN =	ifort
#FORTRAN =	pgf90
#FORTRAN =	pathf90
#FORTRAN =	gfortran
#FORTRAN =	g95
#FORTRAN =	xlf90

F90DOC = ../../bin/f90doc-0.4.0/f90doc

PACKAGE =	CARMA
TGZ =		CARMA.tar

FFLAGS = #-ggdb -O0 -traceback
#FFLAGS += -DSINGLE                    # for single precision
#FFLAGS += -DDEBUG                     # for debug print statements


# Add options for the Intel Fortran compiler.
ifeq ($(FORTRAN),ifort)
#  FFLAGS += -ftz -fp-model precise
  FFLAGS += -fp-model precise
  
  # Work around for an incompatibility with some versions of ifort and OSX.
#  FFLAGS += -use-asm

  # Debug options.
  #FFLAGS += -g -O0 -traceback -fp-stack-check -check bounds -check uninit -fpe0 -ftrapuv
  
  # Open/MP
  # FFLAGS += -qopenmp
  
  # The no_pie flags also the executable to work with idb.
  LDFLAGS = $(FFLAGS) -no_pie
endif

# Add options for the Portland Group compiler.
ifeq ($(FORTRAN),pgf90)
  FFLAGS  += 

  # Debug options.
#  FFLAGS += -g -O0 -Mbounds

  # Open/MP
#  FFLAGS  += -mp

  LDFLAGS = $(FFLAGS)
endif

# Add options for the g95 compiler.
ifeq ($(FORTRAN),g95)
#  FFLAGS  += -fzero -ffree-line-length-huge
  FFLAGS  += -ffree-line-length-huge

  # Debug options.
#  FFLAGS += -g -fbounds-check -ftrace=full
  
  # Open/MP
  #
  # NOTE: g95 does not support Open/MP directives. This will cause one
  # test (carma_test) to fail to link.

  LDFLAGS = $(FFLAGS)
endif

# Add options for the IBM XL Fortran compiler.
#
# NOTE: It doesn't support float to zero.
ifeq ($(FORTRAN),xlf90)
  FFLAGS  += -q64 -qarch=auto -qspillsize=2500 -g -qfullpath 
 
  # Debug options.
  FFLAGS += -qinitauto=7FF7FFFF -qflttrap=ov:zero:inv:en -C
  
  # Open/MP
#  FFLAGS += -qsmp=omp
#  FFLAGS += -qsmp=omp:noopt

  LDFLAGS = $(FFLAGS)
endif

ifeq ($(FORTRAN),gfortran)
  FFLAGS += -fexceptions -fbackslash
  FFLAGS += -fPIC -fno-omit-frame-pointer -ffixed-line-length-none
  #FLIBS= $RPATH $MLIBS -lm
  FOPTIMFLAGS = -O3
  FDEBUGFLAGS = -g

  #LDFLAGS = -pthread -shared -Wl, -Wl,--no-undefined
endif

# Overridning the implicit rules, which would try to use m2c to
# create the .mod.
%.mod : %.o ;
%.o : %.F90 ;
%.html : %.F90 ; $(F90DOC) -cs $<

# Add the directories where the source files are located.
VPATH := ../../source/base ../../tests

# These makefiles have the object lists and dependence information
# for the respective components.
#
# NOTE: In the future it might be nice to generate this dependency
# try automatically.
include ../../source/base/Makefile
include ../../tests/Makefile

# Rules for each executable that could be build.
CDAY.exe : $(CARMA_OBJ) carma_day.o carma_testutils.o atmosphere_mod.o
	$(FORTRAN) $(LDFLAGS) -o CDAY.exe carma_day.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)

CEAST.exe : $(CARMA_OBJ) carma_east.o carma_testutils.o atmosphere_mod.o
	$(FORTRAN) $(LDFLAGS) -o CEAST.exe carma_east.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)

CWEST.exe : $(CARMA_OBJ) carma_west.o carma_testutils.o atmosphere_mod.o
	$(FORTRAN) $(LDFLAGS) -o CWEST.exe carma_west.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)
	
CDAYbig.exe : $(CARMA_OBJ) carma_daybig.o carma_testutils.o atmosphere_mod.o
	$(FORTRAN) $(LDFLAGS) -o CDAYbig.exe carma_daybig.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)

CEASTbig.exe : $(CARMA_OBJ) carma_eastbig.o carma_testutils.o atmosphere_mod.o
	$(FORTRAN) $(LDFLAGS) -o CEASTbig.exe carma_eastbig.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)

CWESTbig.exe : $(CARMA_OBJ) carma_westbig.o carma_testutils.o atmosphere_mod.o
	$(FORTRAN) $(LDFLAGS) -o CWESTbig.exe carma_westbig.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)
	
CDAYsmall.exe : $(CARMA_OBJ) carma_daysmall.o carma_testutils.o atmosphere_mod.o
	$(FORTRAN) $(LDFLAGS) -o CDAYsmall.exe carma_daysmall.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)

CEASTsmall.exe : $(CARMA_OBJ) carma_eastsmall.o carma_testutils.o atmosphere_mod.o
	$(FORTRAN) $(LDFLAGS) -o CEASTsmall.exe carma_eastsmall.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)

CWESTsmall.exe : $(CARMA_OBJ) carma_westsmall.o carma_testutils.o atmosphere_mod.o
	$(FORTRAN) $(LDFLAGS) -o CWESTsmall.exe carma_westsmall.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)

diamondback_test.exe : $(CARMA_OBJ) carma_diamondback_test.o carma_testutils.o atmosphere_mod.o
	$(FORTRAN) $(LDFLAGS) -o diamondback_test.exe carma_diamondback_test.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)


# Compile everything.
all : CDAY.exe CEAST.exe CWEST.exe CDAYbig.exe CEASTbig.exe CWESTbig.exe \
      CDAYsmall.exe CEASTsmall.exe CWESTsmall.exe diamondback_test.exe

# Compile all of the documentation.
doc : $(CARMA_DOC) $(TEST_DOC)

clean:
	/bin/rm -f *.o *.mod *.exe *.txt *.html

# The Mac creates .DS_Store files that we don't want in the tar file, so
# exclude them.
tar:
	tar --directory ../.. -cvf $(TGZ) --exclude .DS_Store --exclude .svn \
	  Makefile make-carma.csh run-carma.csh README run-regress.csh view-bench.csh run-all.csh \
	  source tests bin doc/ChangeLog doc/ChangeLog_template doc/index.html
