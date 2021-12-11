### Convert Precision

This Fortran program converts large binary files in parallel from/to single to/from double precision, using MPI-I/O. The code assumes that the each file contains *only* single- or double-precision floating-point.

### Usage

  1. Build the Fortran program:
```bash
mpif90 -O3 -D_SINGLE_TO_DOUBLE convert_precision.F90 -o single2double # executable for single-to-double conversion
mpif90 -O3 -D_DOUBLE_TO_SINGLE convert_precision.F90 -o double2single # executable for double-to-single conversion
```
n.b.: by default, the compiled executable from `mpif90 convert_precision.F90` will convert from double to single precision.

  2. list the files to be converted in a file `files.in`. For instance:
```bash
ls a.bin b.bin c.bin > files.in
```

  3. run the code, e.g.:
```bash
NUM_TASKS=16
mpirun -n $NUM_TASKS ./double2single
```

  4. done! the files with the same name as the uncoverted ones, are generated with the extension `.converted` appended, i.e., for the example above:
```bash
$ ls *.converted
a.bin.converted b.bin.converted c.bin.converted
```
