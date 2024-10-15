include $(PFUNIT_DIR)/include/PFUNIT.mk

all: test_mods test_core
	./test_mods
	./test_core

%.o : %.F90
	$(FC) -c $(FFLAGS) $<

FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
FFLAGS += -I..

test_core_TESTS := test_bracket.pf test_grid.pf test_adv.pf
test_core_OTHER_LIBRARIES := -L.. -lsut
$(eval $(call make_pfunit_test,test_core))

test_mods_TESTS := test_mod_CO2SYSconstants.pf
test_mods_OTHER_LIBRARIES := -L.. -lsut
$(eval $(call make_pfunit_test,test_mods))

clean:
	$(RM) *.o *.mod *.a  *.inc $(test_mods_TESTS:.pf=.F90) \
		$(test_core_TESTS:.pf=.F90) test_mods test_core

test_core: ../libsut.a
test_mods: ../libsut.a

test_adv.mod: test_grid.mod
