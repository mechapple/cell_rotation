# List of executables
EXECUTABLES= gccm_single gccm_generate

# Makefile rules -L/users/pala2/lore/gmp/lib -lgmp
all: $(EXECUTABLES)

gccm_single: GCCM_SingleCrystal1.f90
	gfortran -o gccm_single GCCM_SingleCrystal1.f90 
	
gccm_generate: GCCM_GenerateCrystal1.f90
	gfortran -o gccm_generate GCCM_GenerateCrystal1.f90

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
