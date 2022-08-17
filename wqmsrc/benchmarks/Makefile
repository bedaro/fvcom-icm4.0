SRCS = co2sys_bench.f90

OBJS = $(SRCS:.f90=.o)

%.o: %.f90
	$(FC) -c $(FFLAGS) $(INCS) -I.. $< -o $@

all: co2sys_bench

co2sys_bench: co2sys_bench.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o co2sys_bench $(CO2SYSDEPS) co2sys_bench.o

CO2SYSDEPS = ../mod_prec.o ../mod_sizes.o ../mod_lims.o ../mod_control.o \
	../mod_hydrovars.o ../mod_wqm.o ../mod_buffers.o ../mod_fileinfo.o \
	../mod_CO2SYSconstants.o

co2sys_bench.o: $(CO2SYSDEPS)

clean:
	rm -f *.o co2sys_bench
