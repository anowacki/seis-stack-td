FC = gfortran
FFLAGS = -g -O2 -fcheck=all -fbacktrace
LDFLAGS = 

# Use OpenMP
FFLAGS += -fopenmp
LDFLAGS += -fopenmp

# Flags for FFTW3
LIB = -L/opt/local/lib -lfftw3f
INC = -I/opt/local/include
FFLAGS += $(INC)
LDFLAGS += $(LIB)

# Flags for NetCDF
LIB += -lnetcdff

DEFS = -DFORCE_BIGENDIAN_SACFILES

B = bin
O = obj
S = src

PROGS = \
   $B/stack_sum \
   $B/stack_vespa

MODS = \
   $O/f90sac.o \
   $O/stack.o

default: progs

progs: $(PROGS)

# Programs
$B/%: $(MODS) $(OBJS) $O/%.o
	$(FC) -o $@ $(LDFLAGS) $^

# Build object files
$O/%.o: $S/%.f90
	$(FC) -o $@ -c $(FFLAGS) -J$O $S/$*.f90

$O/%.o: $S/%.F90
	$(FC) -o $@ $(DEFS) -c $(FFLAGS) -J$O $S/$*.F90

# Make programs
%: $B/%

# Module interdependencies
$O/stack.o: $O/f90sac.o

.SECONDARY:

.PHONY: clean spotless

spotless: clean
	/bin/rm -f $(PROGS)

clean:
	/bin/rm -f $O/*.o $O/*.mod
