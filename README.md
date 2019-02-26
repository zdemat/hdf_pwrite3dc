# hdf_pwrite3dc
writing large series of image-like data into HDF5 dataset

## Requrements
- MPI enabled HDF5
- netCDF
- Szip
- Automake

## Building
```bash
aclocal
automake --add-missing
autoconf

#make clean
#make distclean

LIBS=-lgpfs ./configure --with-mpi=yes IORBSIZE=4M TESTFSPATH=/tmp --enable-shared HDF5_USE_SHLIB=yes
make modules
make HDF5_USE_SHLIB=yes
```
