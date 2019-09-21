This is the code corresponding to the paper:

> E. G. Birgin, J. L. Gardenghi, J. M. Martínez, and S. A. Santos, On
  the use of third-order models with fourth-order regularization for
  unconstrained optimization, 2018.

This distribution includes Algencan 3.1.1.

## Compiling and running the code

Consider $(AR4) the folder this README is located.

1. Go to $(AR4) folder and type
```
make
```

2. You can run the drivers to MGH problems and the Packing problem in
`$(AR4)/bin` folder by running `mgh` and `pack` binaries files,
respectively. Follow the instructions given by these programs to run
the tests.

## Algencan and HSL

The tests in the paper were run using MA57 and BLAS
subroutines. Algencan does not need these files to work, but in order
to reproduce the tests in the paper, the user must include them. To
include them, follow these steps:

1. Obtain HSL_MA57 from Harwell Subroutine Library
[website](http://www.hsl.rl.ac.uk/).

2. Obtain BLAS subroutine from Netlib
[website](http://www.netlib.org/blas/).

3. Copy the files in $(AR4)/algencan/sources/hsl folder:

   * From BLAS library: `dgemm.f`, `dgemv.f`, `dtpmv.f`, `dtpsv.f`,
     `idamax.f', `lsame.f`, `xerbla.f`.

   * From HSL MA57 library: `hsl_ma57d.f90`, `hsl_zd11d.f90`,
     `ma57ad.f`, `mc21ad.f`, `mc34ad.f`, `mc47ad.f`, `mc59ad.f`,
     `mc64ad.f`, `mc71ad.f`. (Except `hsl_ma57d.f90`, all other
     routines usually comes inside files named `deps.f` and
     `deps90.f90`. The user may need to save them in files with these
     names.)

   * HSL MA57 requires MeTIS. You can use `fakemetis.f` that comes
     with HSL MA57 distribution.

4. Include `-lhsl` flag in variable `ALGENCAN_FLAGS` inside
`$(AR4)/Makefile` file.