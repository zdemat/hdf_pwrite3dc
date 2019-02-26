AUTOMAKE_OPTIONS = foreign
SUBDIRS = src

local_ior = ior-3.0.1

bin_PROGRAMS = hdf_pwrite3dc
EXTRA_DIST = conf tests results doc $(local_ior).tar.gz
hdf_pwrite3dc_SOURCES = src/hdf_pwrite3dc.f90 src/kinds.f90 $(local_ior).tar.gz

foo.$(OBJEXT): src/kinds.f90

#kinds.mod: kinds

#src/kinds.mod: src/kinds.f90
#	cd src; $(FC) $(FFLAGS) $(LDFLAGS) -c kinds.f90

mpirun:
	echo "Creating symbolic link to mpirun."
	ln -s $(MPICMD) mpirun

clean-local:
	-rm -rf conf ior mpirun
	-rm -rfd $(local_ior)
	-rm -rf src/*.mod

conf:
	mkdir conf

results:
	mkdir results

ior:
	echo "Building ior."
	rm -rfd $(local_ior)
	tar -xf $(srcdir)/$(local_ior).tar.gz
	cd $(local_ior); pwd; ./bootstrap
	cd $(local_ior); ./configure --with-mpiio --with-hdf5 LIBS=$(LIBS) CC=h5pcc CFLAGS=-DH5_USE_16_API
	cd $(local_ior); make
	echo "Creating symbolic link to ior."
	ln -s $(local_ior)/src/ior ior

all: hdf_pwrite3dc ior conf mpirun

#FORCE: # no dependencies or commands
#.PHONY: FORCE # hint for gnu make, ignored by others
