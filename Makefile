# Compiler and flags
FC = gfortran
FLAGS = -cpp -Wall -O3 -DOFFLINE

# NetCDF-Fortran installation. Default: serial gcc build on DKRZ Levante.
# Override on the command line, e.g. make NETCDF_FORTRAN=/path/to/netcdf-fortran
NETCDF_FORTRAN ?= /sw/spack-levante/netcdf-fortran-4.5.3-l2ulgp

# NetCDF library and include paths (rpath so no LD_LIBRARY_PATH is needed at runtime)
NETCDF_LIB ?= $(shell $(NETCDF_FORTRAN)/bin/nf-config --flibs) -Wl,-rpath,$(NETCDF_FORTRAN)/lib
NETCDF_INC ?= -I$(NETCDF_FORTRAN)/include
# Source files
SRC = dust_only.f90 src_dust.f90 data_dust.f90 mo_dust.f90

# Object and module directory
OBJDIR = objs

# Object files
OBJ = $(addprefix $(OBJDIR)/,$(SRC:.f90=.o))
MOD = $(addprefix $(OBJDIR)/,$(SRC:.f90=.mod))

# Executable
EXECUTABLE = dustonly

# Main target
all: $(EXECUTABLE)

# Rule to build the executable
$(EXECUTABLE): $(OBJ)
	$(FC) $(FLAGS) -o $@ $^ $(NETCDF_LIB) $(LDFLAGS)

# Rules to compile the source files
$(OBJDIR)/%.o: %.f90
	@mkdir -p $(OBJDIR)
	$(FC) $(FLAGS) -c -o $@ $< $(NETCDF_INC) -I ./$(OBJDIR) $(CPPFLAGS)
	@if [ `ls -1 *.mod 2>/dev/null | wc -l` != 0 ]; then mv *.mod $(OBJDIR)/; fi
	

# Dependencies
$(OBJDIR)/dust_only.o: $(OBJDIR)/src_dust.o $(OBJDIR)/data_dust.o $(OBJDIR)/mo_dust.o
$(OBJDIR)/src_dust.o: $(OBJDIR)/data_dust.o $(OBJDIR)/mo_dust.o
$(OBJDIR)/data_dust.o: $(OBJDIR)/mo_dust.o

# Clean target
clean:
	rm -f $(EXECUTABLE) $(OBJ) $(MOD)
	rm -r $(OBJDIR)

