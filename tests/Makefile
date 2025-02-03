# Sub makefile for test files

# This is intended to be included by another makefile to actually
# build the system. It has all the dependency information for the
# files in the test tree.

# NOTE: In the future (or in parnet models) this could be replaced
# by automatic dependency generation.

TEST_OBJ = carma_day.o carma_east.o carma_west.o carma_daybig.o carma_eastbig.o \
           carma_westbig.o carma_daysmall.o carma_eastsmall.o carma_westsmall.o
           

TEST_DOC = carma_day.html carma_east.html carma_west.html carma_daybig.html carma_eastbig.html \
           carma_westbig.html carma_daysmall.html carma_eastsmall.html carma_westsmall.html

atmosphere_mod.o : atmosphere_mod.F90
	$(FORTRAN) $(FFLAGS) -c $<

carma_testutils.o : carma_testutils.F90
	$(FORTRAN) $(FFLAGS) -c $<

carma_day.o : carma_day.F90 atmosphere_mod.mod carma_mod.mod carma_planet_mod.mod carma_condensate_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

carma_east.o : carma_east.F90 atmosphere_mod.mod carma_mod.mod carma_planet_mod.mod carma_condensate_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

carma_west.o : carma_west.F90 atmosphere_mod.mod carma_mod.mod carma_planet_mod.mod carma_condensate_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

carma_daybig.o : carma_daybig.F90 atmosphere_mod.mod carma_mod.mod carma_planet_mod.mod carma_condensate_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

carma_eastbig.o : carma_eastbig.F90 atmosphere_mod.mod carma_mod.mod carma_planet_mod.mod carma_condensate_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

carma_westbig.o : carma_westbig.F90 atmosphere_mod.mod carma_mod.mod carma_planet_mod.mod carma_condensate_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<
	
carma_daysmall.o : carma_daysmall.F90 atmosphere_mod.mod carma_mod.mod carma_planet_mod.mod carma_condensate_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

carma_eastsmall.o : carma_eastsmall.F90 atmosphere_mod.mod carma_mod.mod carma_planet_mod.mod carma_condensate_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

carma_westsmall.o : carma_westsmall.F90 atmosphere_mod.mod carma_mod.mod carma_planet_mod.mod carma_condensate_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

carma_diamondback_test.o : carma_diamondback_test.F90 atmosphere_mod.mod carma_mod.mod carma_planet_mod.mod carma_condensate_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<