include $(PFUNIT_DIR)/include/PFUNIT.mk

all: test_mods
	./test_mods

%.o : %.F90
	$(FC) -c $(FFLAGS) $<

FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
FFLAGS += -I..

test_mods_TESTS := test_mod_CO2SYSconstants.pf
test_mods_OTHER_LIBRARIES := -L.. -lsut
$(eval $(call make_pfunit_test,test_mods))

clean:
	$(RM) *.o *.mod *.a  *.inc $(test_mods_TESTS:.pf=.F90)
